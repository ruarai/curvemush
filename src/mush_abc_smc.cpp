
#include "musher.h"
#include "mush_params.h"
#include "rmultinom_vec.h"

#include "abc_smc.h"
#include "abc_smc_data.h"

#include <random>

#include <RcppArmadillo.h>
#include <RcppThread.h>


using namespace Rcpp;


// [[Rcpp::export]]
List mush_abc_smc(
    int n_parameter_samples,
    int n_particles,
    int n_delay_samples,

    int n_days,
    int steps_per_day,

    NumericMatrix ensemble_curves,

    NumericMatrix mat_pr_age_given_case,
    NumericMatrix mat_pr_hosp,
    NumericMatrix mat_pr_ICU,

    DataFrame forecasting_parameters,

    std::vector<float> thresholds_vec,
    
    std::vector<int> known_ward_vec,
    std::vector<int> known_ICU_vec,

    int day_start_fit,
    int day_end_fit,

    arma::mat spline_basis
) {
    int n_strat_samples = forecasting_parameters.nrow() / def_n_strats;

    smc_inputs input_data {
        .case_curves = std::vector<std::vector<std::vector<int>>>(def_n_strats, std::vector<std::vector<int>>(n_parameter_samples, std::vector<int>(n_days, 0))),
        .pr_hosp_curves = std::vector<std::vector<std::vector<float>>>(def_n_strats, std::vector<std::vector<float>>(n_parameter_samples, std::vector<float>(n_days, -1))),
        .pr_ICU_curves = std::vector<std::vector<std::vector<float>>>(def_n_strats, std::vector<std::vector<float>>(n_parameter_samples, std::vector<float>(n_days, -1))),

        .strat_datas = std::vector<std::vector<strat_data>>(def_n_strats, std::vector<strat_data>(n_strat_samples)),
    };
    
    input_data.n_days = n_days;
    input_data.steps_per_day = steps_per_day;
    input_data.n_delay_samples = n_delay_samples;

    input_data.n_particles = n_particles;

    input_data.n_parameter_samples = n_parameter_samples;
    input_data.n_strat_samples = n_strat_samples;

    input_data.day_start_fit = day_start_fit;
    input_data.day_end_fit = day_end_fit;
    
    input_data.known_ward_occupancy = known_ward_vec;
    input_data.known_ICU_occupancy = known_ICU_vec;

    input_data.thresholds = thresholds_vec;

    input_data.spline_basis = spline_basis;

    for(int s = 0; s < def_n_strats; s++) {
        for(int i = 0; i < n_strat_samples; i++)
            input_data.strat_datas[s][i] = strat_data::read_strat_data(forecasting_parameters, i * def_n_strats + s);
    }

    int n_curves = ensemble_curves.ncol();
    // Producing n_parameter_samples case/hospitalisation trajectories to feed into the model:
    for (int i = 0; i < n_parameter_samples; i++)
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
                input_data.case_curves[s][i][d] = age_samples[s];

            
        }

        for(int d = 0; d < n_days; d++) {
            for (int s = 0; s < def_n_strats; s++) {
                input_data.pr_hosp_curves[s][i][d] = mat_pr_hosp(d * def_n_strats + s, i % mat_pr_hosp.ncol());
                input_data.pr_ICU_curves[s][i][d] = mat_pr_ICU(d * def_n_strats + s, i % mat_pr_ICU.ncol());
            }
        }
    }

    smc_results smc_out = abc_smc::process_abc_smc(
        input_data
    );

    int result_size = input_data.n_days * def_n_compartments;
    int n_results = n_particles * result_size;

    NumericVector sample_label(n_results);
    NumericVector t_day(n_results);
    NumericVector compartment_label(n_results);
    NumericVector compartment_count(n_results);
    NumericVector compartment_transitions(n_results);

    for (int i = 0; i < n_particles; i++)
    {
        for (int x = 0; x < result_size; x++)
        {
            int ix = i * result_size + x;

            sample_label[ix] = i;

            t_day[ix] = x % n_days;
            compartment_label[ix] = smc_out.results[i][0].occupancy_compartment_labels[x];

            compartment_count[ix] = 0;
            compartment_transitions[ix] = 0;
            for (int s = 0; s < def_n_strats; s++)
            {
                compartment_count[ix] += smc_out.results[i][s].occupancy_counts[x];
                compartment_transitions[ix] += smc_out.results[i][s].transitions[x];
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
    int grped_result_size = n_days * def_n_compartment_groups;
    int grped_n_results = n_particles * grped_result_size;

    NumericVector grped_sample_label(grped_n_results);
    NumericVector grped_t_day(grped_n_results);
    NumericVector grped_compartment_group_label(grped_n_results);
    NumericVector grped_compartment_group_count(grped_n_results);
    NumericVector grped_compartment_group_transitions(grped_n_results);



    for (int i = 0; i < n_particles; i++)
    {
        for (int x = 0; x < grped_result_size; x++)
        {
            int ix = i * grped_result_size + x;

            grped_sample_label[ix] = i;
            grped_t_day[ix] = x % n_days;
            grped_compartment_group_label[ix] = smc_out.results[i][0].grouped_occupancy_compartment_labels[x];

            grped_compartment_group_count[ix] = 0;
            grped_compartment_group_transitions[ix] = 0;
            for (int s = 0; s < def_n_strats; s++)
            {
                grped_compartment_group_count[ix] += smc_out.results[i][s].grouped_occupancy_counts[x];
                grped_compartment_group_transitions[ix] += smc_out.results[i][s].grouped_transitions[x];
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

    return List::create(
        _["results"] = R_results,
        _["grouped_results"] = R_grped_results,

        _["smc_weights"] = smc_out.weights,
        _["smc_attempts"] = smc_out.attempts,

        _["smc_spline_fits"] = smc_out.spline_fits
    );
}