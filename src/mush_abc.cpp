

#include <random>

#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"
#include "rmultinom_vec.h"

using namespace Rcpp;

float logit(float x) { return std::log(x / (1 - x)); }
float ilogit(float x) { return std::exp(x) / (std::exp(x) + 1);}

// [[Rcpp::export]]
List mush_abc(
    int n_samples,
    int n_delay_samples,

    int n_outputs,

    int n_days,
    int steps_per_day,

    int t_forecast_start,

    int ward_threshold,
    int ICU_threshold,

    NumericMatrix ensemble_curves,

    NumericMatrix mat_pr_age_given_case,
    NumericMatrix mat_pr_hosp,
    NumericMatrix mat_pr_ICU,

    DataFrame forecasting_parameters,
    
    IntegerVector known_ward_vec,
    IntegerVector known_ICU_vec,
    
    double prior_sigma_los,
    double prior_sigma_hosp)
{
    mush_params params;
    params.n_days = n_days;
    params.steps_per_day = steps_per_day;
    params.n_delay_samples = n_delay_samples;

    int n_strat_samples = forecasting_parameters.nrow() / def_n_strats;
    std::vector<strat_data> strat_datas(n_samples * def_n_strats);

    for (int i = 0; i < n_samples * def_n_strats; i++)
        strat_datas[i] = strat_data::read_strat_data(forecasting_parameters, i);



    // Define and initialize our hospitalisation curves (vector of vectors of vectors)
    std::vector<std::vector<std::vector<int>>> case_curves(
        def_n_strats,
        std::vector<std::vector<int>>(n_samples, std::vector<int>(n_days, 0)));

    std::vector<std::vector<std::vector<float>>> pr_hosp_curves(
        def_n_strats,
    std::vector<std::vector<float>>(n_samples, std::vector<float>(n_days, -1)));

    std::vector<std::vector<std::vector<float>>> pr_ICU_curves(
        def_n_strats,
        std::vector<std::vector<float>>(n_samples, std::vector<float>(n_days, -1)));

    int n_curves = ensemble_curves.ncol();
    // Producing n_samples case/hospitalisation trajectories to feed into the model:
    for (int i = 0; i < n_samples; i++)
    {
        // Select a curve from our backcast/ensemble curve matrix
        // Looping (by % n_curves) if we are performing more samples than we have ensembles
        NumericVector curve_i = ensemble_curves(_, i % n_curves);

        // Over the backcast period, cases are split by age group (and known precisely)
        for (int s = 0; s < def_n_strats; s++)
            for (int d = 0; d < t_forecast_start; d++)
                case_curves[s][i][d] = curve_i[d * def_n_strats + s];

        // But in the forecast period, we perform the age sampling ourselves
        for (int d = 0; d < n_days - t_forecast_start; d++)
        {
            int n_cases = curve_i[t_forecast_start * def_n_strats + d];

            if (n_cases == 0)
                continue;

            NumericVector pr_age_given_case(def_n_strats);
            for(int s = 0; s < def_n_strats; s++)
                pr_age_given_case[s] = mat_pr_age_given_case(d * def_n_strats + s, i % mat_pr_age_given_case.ncol());

            // Produce a vector distributing the n_cases across each age stratification
            // according to pr_age_given_case
            IntegerVector age_samples = rmultinom_vec(n_cases, pr_age_given_case);

            for (int s = 0; s < def_n_strats; s++)
                case_curves[s][i][t_forecast_start + d] = age_samples[s];

            
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
        

    RcppThread::parallelFor(
        0, n_outputs,
        [&results, params, &case_curves, strat_datas, n_strat_samples, known_ward, known_ICU, n_samples, t_forecast_start,
         ward_threshold, ICU_threshold, &los_scale_samples, &pr_hosp_scale_samples, &prior_chosen,
         &pr_hosp_curves, &pr_ICU_curves](unsigned int i)
        {
            bool rejected = true;
            std::random_device rd;  
            std::mt19937 rng(rd()); 
            std::uniform_int_distribution<int> uni(0, n_samples - 1); 

            while(rejected) {

                std::vector<mush_results> sample_results(def_n_strats);

                int i_prior = uni(rng);

                        
                for (int s = 0; s < def_n_strats; s++)
                {
                    strat_data s_data = strat_datas[i_prior % n_strat_samples + s];


                    std::vector<int> case_curve_i = case_curves[s][i_prior];

                    for (int d = 0; d < params.n_days - t_forecast_start; d++)
                    {
                        int n_cases = case_curve_i[t_forecast_start + d];

                        if (n_cases == 0)
                            continue;

                        float pr_hosp = pr_hosp_curves[s][i_prior][d];

                        float effect_adj = 1 - std::min(t_forecast_start - d, 60) / 60;

                        float adj_pr_hosp = ilogit(logit(pr_hosp) + pr_hosp_scale_samples[i_prior] * effect_adj);
                            
                        std::binomial_distribution<int> hosp_binom(n_cases, adj_pr_hosp);

                        case_curve_i[t_forecast_start + d] = hosp_binom(rng);
                    }


                    
                    sample_results[s] = musher::mush_curve(
                        params,
                        case_curve_i,
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
                    if(known_ward[t] != -1 && std::abs(ward_counts[t] - known_ward[t]) > ward_threshold) {
                        rejected = true;
                    }
                    if(known_ICU[t] != -1 && std::abs(ICU_counts[t] - known_ICU[t]) > ICU_threshold) {
                        rejected = true;
                    }
                }

                if(!rejected) {
                    results[i] = sample_results;
                    prior_chosen[i] = i_prior;
                }
            }
        }
    );

    int result_size = params.n_days * def_n_compartments;
    int n_results = n_outputs * result_size;

    NumericVector sample_label(n_results);
    NumericVector t_day(n_results);
    NumericVector compartment_label(n_results);
    NumericVector compartment_count(n_results);
    NumericVector compartment_transitions(n_results);
    NumericVector sample_los_scale(n_results);

    for (int i = 0; i < n_outputs; i++)
    {
        for (int x = 0; x < result_size; x++)
        {
            int ix = i * result_size + x;

            sample_label[ix] = i;
            sample_los_scale[ix] = los_scale_samples[i];

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
        _["transitions"] = compartment_transitions,
        _["los_scale"] = sample_los_scale
    );

    // Build out our result dataframe to return
    int grped_result_size = params.n_days * def_n_compartment_groups;
    int grped_n_results = n_outputs * grped_result_size;

    NumericVector grped_sample_label(grped_n_results);
    NumericVector grped_t_day(grped_n_results);
    NumericVector grped_compartment_group_label(grped_n_results);
    NumericVector grped_compartment_group_count(grped_n_results);
    NumericVector grped_compartment_group_transitions(grped_n_results);
    NumericVector grped_sample_los_scale(grped_n_results);
    NumericVector grped_sample_pr_hosp_scale(grped_n_results);


    for (int i = 0; i < n_outputs; i++)
    {
        for (int x = 0; x < grped_result_size; x++)
        {
            int ix = i * grped_result_size + x;

            grped_sample_label[ix] = i;
            grped_t_day[ix] = x % params.n_days;
            grped_compartment_group_label[ix] = results[i][0].grouped_occupancy_compartment_labels[x];
            
            grped_sample_los_scale[ix] = los_scale_samples[prior_chosen[i]];
            grped_sample_pr_hosp_scale[ix] = pr_hosp_scale_samples[prior_chosen[i]];

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
        _["transitions"] = grped_compartment_group_transitions,
        _["los_scale"] = grped_sample_los_scale,
        _["pr_hosp_scale"] = grped_sample_pr_hosp_scale
    );

    List R_list_results = List::create(
        _["results"] = R_results,
        _["grouped_results"] = R_grped_results);

    return R_list_results;
}