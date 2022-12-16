

#include <random>

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"
#include "rmultinom_vec.h"

#include "hosp_sampling.h"

using namespace Rcpp;

//' Perform simulation and optionally fitting
//' @export
// [[Rcpp::export]]
List mush_abc(
    int n_samples,
    int n_delay_samples,

    int n_outputs,

    int n_days,
    int steps_per_day,

    NumericVector thresholds_vec,
    int rejections_per_selections,
    bool do_ABC,

    NumericMatrix ensemble_curves,

    NumericMatrix mat_pr_age_given_case,
    NumericMatrix mat_pr_hosp,
    NumericMatrix mat_pr_ICU,

    DataFrame forecasting_parameters,
    
    IntegerVector known_ward_vec,
    IntegerVector known_ICU_vec,
    
    double prior_sigma_los,
    double prior_sigma_hosp
) {
    mush_params params;
    params.n_days = n_days;
    params.steps_per_day = steps_per_day;
    params.n_delay_samples = n_delay_samples;

    std::vector<strat_data> strat_datas(n_samples * def_n_strats);

    for (int i = 0; i < n_samples * def_n_strats; i++)
        strat_datas[i] = strat_data::read_strat_data(forecasting_parameters, i % forecasting_parameters.nrow());


    // Define and initialize our hospitalisation curves (vector of vectors of vectors)
    std::vector<std::vector<std::vector<int>>> case_curves(def_n_strats, std::vector<std::vector<int>>(n_samples, std::vector<int>(n_days, 0)));

    // pr_hosp and pr_ICU bootstrap curves (to be filled from mat_pr_hosp and mat_pr_ICU)
    std::vector<std::vector<std::vector<float>>> pr_hosp_curves(def_n_strats, std::vector<std::vector<float>>(n_samples, std::vector<float>(n_days, -1)));
    std::vector<std::vector<std::vector<float>>> pr_ICU_curves(def_n_strats, std::vector<std::vector<float>>(n_samples, std::vector<float>(n_days, -1)));

    int n_curves = ensemble_curves.ncol();
    // Producing n_samples case/hospitalisation trajectories to feed into the model:
    for (int i = 0; i < n_samples; i++)
    {
        // Select a curve from our backcast/ensemble curve matrix
        // Looping (by % n_curves) if we are performing more samples than we have ensembles
        NumericVector curve_i = ensemble_curves(_, i % n_curves);

        // Perform the age sampling ourselves
        for (int d = 0; d < n_days; d++)
        {
            int n_cases = curve_i[d];

            if (n_cases == 0)
                continue;

            bool all_zero_pr = true;

            NumericVector pr_age_given_case(def_n_strats);
            for(int s = 0; s < def_n_strats; s++) {
                pr_age_given_case[s] = mat_pr_age_given_case(d * def_n_strats + s, i % mat_pr_age_given_case.ncol());

                if(pr_age_given_case[s] > 0.0001)
                    all_zero_pr = false;
            }

            // There's no probability of cases by any age, skip (otherwise the multinomial sampling will be unhappy)
            if(all_zero_pr)
                continue;

            // Produce a vector distributing the n_cases across each age stratification
            // according to pr_age_given_case
            IntegerVector age_samples = rmultinom_vec(n_cases, pr_age_given_case);

            for (int s = 0; s < def_n_strats; s++)
                case_curves[s][i][d] = age_samples[s];

            
        }

        for(int d = 0; d < n_days; d++) {
            for (int s = 0; s < def_n_strats; s++) {
                pr_hosp_curves[s][i][d] = mat_pr_hosp(d * def_n_strats + s, i % mat_pr_hosp.ncol());
                pr_ICU_curves[s][i][d] = mat_pr_ICU(d * def_n_strats + s, i % mat_pr_ICU.ncol());
            }
        }
    }

    std::vector<float> los_scale_samples(n_samples);
    std::vector<float> pr_hosp_scale_samples(n_samples);
    std::vector<int> known_ward = as<std::vector<int>>(known_ward_vec);
    std::vector<int> known_ICU = as<std::vector<int>>(known_ICU_vec);

    for (int i = 0; i < los_scale_samples.size(); i++)
        los_scale_samples[i] = R::rnorm(0, prior_sigma_los);
        
    for (int i = 0; i < pr_hosp_scale_samples.size(); i++)
        pr_hosp_scale_samples[i] = R::rnorm(0, prior_sigma_hosp);

    std::vector<std::vector<mush_results>> results(n_outputs);
    std::vector<int> prior_chosen(n_outputs);
        
    std::vector<float> thresholds = as<std::vector<float>>(thresholds_vec);
    std::vector<int> n_rejected(thresholds.size(), 0);
    std::vector<int> n_accepted(thresholds.size(), 0);
    int max_rejections = n_outputs * rejections_per_selections;

    // Iterate over the thresholds in thresholds_vec
    for(int i_threshold = 0; i_threshold < thresholds.size(); i_threshold++) {
        

        // For each threshold, try and produce 'n_outputs' results, with rejection of trajectories
        // that do not match true ICU or ward occupancy counts +- the threshold (proportionally)
        RcppThread::parallelFor(
            0, n_outputs,
            [&results, params, &case_curves, strat_datas, known_ward, known_ICU, n_samples,
            &los_scale_samples, &pr_hosp_scale_samples, &prior_chosen, &n_rejected, &n_accepted, thresholds,
             i_threshold, max_rejections, do_ABC, &pr_hosp_curves, &pr_ICU_curves](unsigned int i)
            {

                bool rejected = true;
                std::random_device rd;  
                std::mt19937 rng(rd()); 
                std::uniform_int_distribution<int> uni(0, n_samples - 1);

                // This will repeat until a good trajectory (and corresponding prior) is found
                // If this process -- across all threads -- fails max_rejections times, we will move to the next bigger threshold
                while(rejected) {
                    if(n_rejected[i_threshold] >= max_rejections)
                        return;

                    std::vector<mush_results> sample_results(def_n_strats);
                    int i_prior = uni(rng);

                    for (int s = 0; s < def_n_strats; s++)
                    {

                        strat_data s_data = strat_datas[i_prior * def_n_strats + s];


                        std::vector<int> hospitalised_cases = sample_hospitalised_cases(
                            case_curves[s][i_prior],
                            pr_hosp_curves[s][i_prior],

                            pr_hosp_scale_samples[i_prior],

                            rng
                        );
                        
                        sample_results[s] = musher::mush_curve(
                            params,
                            hospitalised_cases,
                            s_data,
                            los_scale_samples[i_prior],
                            pr_ICU_curves[s][i_prior]
                        );
                    }
                    
                    std::vector<int> ward_counts(params.n_days, 0);
                    std::vector<int> ICU_counts(params.n_days, 0);

                    for(int s = 0; s < def_n_strats; s++) {
                        for(int x = 0; x < params.n_days * def_n_compartment_groups; x++) {

                            if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ward) {
                                ward_counts[x % params.n_days] += sample_results[s].grouped_occupancy_counts[x];
                            } else if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ICU) {
                                ICU_counts[x % params.n_days] += sample_results[s].grouped_occupancy_counts[x];
                            }
                        }
                        
                    }

                    
                    rejected = false;

                    for(int t = 0; t < params.n_days; t++) {
                        if(known_ward[t] != -1 && 
                            std::abs(ward_counts[t] - known_ward[t]) > std::max(known_ward[t] * thresholds[i_threshold], 2.0f)) {
                            rejected = true;
                        }
                        if(known_ICU[t] != -1 && 
                            std::abs(ICU_counts[t] - known_ICU[t]) > std::max(known_ICU[t] * thresholds[i_threshold] * 1.5f, 4.0f)) {
                            rejected = true;
                        }
                    }

                    if(!rejected || !do_ABC) {
                        results[i] = sample_results;
                        prior_chosen[i] = i_prior;
                        n_accepted[i_threshold]++;

                        rejected = false;
                    } else {
                        n_rejected[i_threshold]++;
                    }
                }
            }
        );

        if(n_rejected[i_threshold] < max_rejections || !do_ABC)
            break;
    }




    int result_size = params.n_days * def_n_compartments;
    int n_results = n_outputs * result_size;

    NumericVector sample_label(n_results);
    NumericVector t_day(n_results);
    NumericVector compartment_label(n_results);
    NumericVector compartment_count(n_results);
    NumericVector compartment_transitions(n_results);

    for (int i = 0; i < n_outputs; i++)
    {
        for (int x = 0; x < result_size; x++)
        {
            int ix = i * result_size + x;

            sample_label[ix] = i;

            t_day[ix] = x % params.n_days;
            compartment_label[ix] = results[i][0].occupancy_compartment_labels[x];

            compartment_count[ix] = 0;
            compartment_transitions[ix] = 0;
            for (int s = 0; s < def_n_strats; s++)
            {
                compartment_count[ix] += results[i][s].occupancy_counts[x];
                compartment_transitions[ix] += results[i][s].transitions[x];
            }
        }
    }

    DataFrame R_results = DataFrame::create(
        _["sample"] = sample_label,
        _["t_day"] = t_day,
        _["compartment"] = compartment_label,
        _["count"] = compartment_count,
        _["transitions"] = compartment_transitions
    );

    // Build out our result dataframe to return
    int grped_result_size = params.n_days * def_n_compartment_groups;
    int grped_n_results = n_outputs * grped_result_size;

    NumericVector grped_sample_label(grped_n_results);
    NumericVector grped_t_day(grped_n_results);
    NumericVector grped_compartment_group_label(grped_n_results);
    NumericVector grped_compartment_group_count(grped_n_results);
    NumericVector grped_compartment_group_transitions(grped_n_results);

    for (int i = 0; i < n_outputs; i++)
    {
        for (int x = 0; x < grped_result_size; x++)
        {
            int ix = i * grped_result_size + x;

            grped_sample_label[ix] = i;
            grped_t_day[ix] = x % params.n_days;
            grped_compartment_group_label[ix] = results[i][0].grouped_occupancy_compartment_labels[x];

            grped_compartment_group_count[ix] = 0;
            grped_compartment_group_transitions[ix] = 0;
            for (int s = 0; s < def_n_strats; s++)
            {
                grped_compartment_group_count[ix] += results[i][s].grouped_occupancy_counts[x];
                grped_compartment_group_transitions[ix] += results[i][s].grouped_transitions[x];
            }
        }
    }

    DataFrame R_grped_results = DataFrame::create(
        _["sample"] = grped_sample_label,
        _["t_day"] = grped_t_day,
        _["compartment_group"] = grped_compartment_group_label,
        _["count"] = grped_compartment_group_count,
        _["transitions"] = grped_compartment_group_transitions
    );

    
    // Build out our result dataframe to return
    int aged_grped_n_results = n_outputs * grped_result_size * def_n_strats;

    NumericVector aged_grped_sample_label(aged_grped_n_results);
    NumericVector aged_grped_t_day(aged_grped_n_results);
    NumericVector aged_grped_age_group(aged_grped_n_results);
    NumericVector aged_grped_compartment_group_label(aged_grped_n_results);
    NumericVector aged_grped_compartment_group_count(aged_grped_n_results);
    NumericVector aged_grped_compartment_group_transitions(aged_grped_n_results);

    int ix = 0;
    for (int i = 0; i < n_outputs; i++)
    {
        for (int x = 0; x < grped_result_size; x++)
        {
            for (int s = 0; s < def_n_strats; s++)
            {
                aged_grped_sample_label[ix] = i;
                aged_grped_t_day[ix] = x % params.n_days;
                aged_grped_compartment_group_label[ix] = results[i][s].grouped_occupancy_compartment_labels[x];
                aged_grped_age_group[ix] = s;

                aged_grped_compartment_group_count[ix] = results[i][s].grouped_occupancy_counts[x];
                aged_grped_compartment_group_transitions[ix] = results[i][s].grouped_transitions[x];

                ix++;
            }
        }
    }

    DataFrame R_aged_grped_results = DataFrame::create(
        _["sample"] = aged_grped_sample_label,
        _["t_day"] = aged_grped_t_day,
        _["compartment_group"] = aged_grped_compartment_group_label,
        _["age_group"] = aged_grped_age_group,
        _["count"] = aged_grped_compartment_group_count,
        _["transitions"] = aged_grped_compartment_group_transitions
    );

    return List::create(
        _["results"] = R_results,
        _["grouped_results"] = R_grped_results,
        _["age_stratified_grouped_results"] = R_aged_grped_results,

        _["n_rejected"] = n_rejected,
        _["n_accepted"] = n_accepted,

        _["prior_chosen"] = prior_chosen,
        _["prior_los_scale"] = los_scale_samples,
        _["prior_pr_hosp"] = pr_hosp_scale_samples
    );
}