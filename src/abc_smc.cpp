#include "abc_smc.h"
#include "abc_smc_data.h"

#include "hosp_sampling.h"

#include <RcppArmadillo.h>
#include <RcppThread.h>
using namespace Rcpp;

#define p_pr_hosp_offset 0
#define p_los_offset 1

typedef std::vector<double> stdvec;


double normal_pdf(double x, double m, double s) {
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

float vec_normalised_rmse(arma::vec &x, arma::vec &x_est) {
    arma::vec square_diff = arma::square(x - x_est);

    float rmse = std::sqrt(arma::sum(square_diff) / x.size());
    float mean = arma::mean(x);

    return rmse / mean;
}

smc_results abc_smc::process_abc_smc(
    smc_inputs &input_data
) {

    int n_particles = input_data.n_particles;
    int n_thresholds = input_data.thresholds.size();

    int n_params = 2;

    std::vector<arma::mat> params(n_particles, arma::mat(n_params, n_thresholds, arma::fill::zeros));

    arma::mat weights(n_particles, n_thresholds);
    arma::mat attempts(n_particles, n_thresholds);


    int n_fit_days = 0;
    int day_start_fit = -1;
    for(int d = 0; d < input_data.n_days; d++) {
        if(input_data.known_ward_occupancy[d] != -1) {
            n_fit_days++;

            if(day_start_fit == -1)
                day_start_fit = d;
        }
    }
    
    arma::vec known_ward_occupancy(n_fit_days);
    arma::vec known_ICU_occupancy(n_fit_days);
    for(int d = 0; d < n_fit_days; d++) {
        known_ward_occupancy(d) = input_data.known_ward_occupancy[d + day_start_fit];
        known_ICU_occupancy(d) = input_data.known_ICU_occupancy[d + day_start_fit];
    }

    
    std::vector<std::vector<mush_results>> results(n_particles);

    double perturb_sigma = 0.025;

    for(int t = 0; t < input_data.thresholds.size(); t++) {
		Rcout << "t = " << t << ", (" << input_data.thresholds[t] << ")\n";

        RcppThread::parallelFor(0, n_particles,
        
        [&params, &weights, &attempts, &results, &input_data, &known_ward_occupancy, &known_ICU_occupancy,
         n_params, n_fit_days, day_start_fit, perturb_sigma, n_particles, t] (unsigned int p) {
            std::random_device rd;
			std::mt19937 gen(rd());

			std::normal_distribution<> perturb_kernel(0, perturb_sigma);
			std::normal_distribution<> prior(0, 1);

            std::uniform_int_distribution<int> param_sample_dist(0, input_data.n_parameter_samples - 1);


            bool fit_achieved = false;
            int n_attempts = 0;

            std::vector<mush_results> sample_results;
            arma::vec spline_params;
            int i_prior = -1;

            while(!fit_achieved) {
                
				if(t == 0) { // Sample from prior at t == 0
					for (int w = 0; w < n_params; w++)
						params[p](w, t) = prior(gen);
				}
				else {
					std::vector<double> prev_weights = arma::conv_to<stdvec>::from(weights.col(t - 1));
					std::discrete_distribution<> sample_dist(prev_weights.begin(), prev_weights.end());

					int p_sample = sample_dist(gen);

					for(int w = 0; w < n_params; w++)
						params[p](w, t) = params[p_sample](w, t - 1) + perturb_kernel(gen);
					
				}

                spline_params = params[p].col(t);
                
                
                i_prior = param_sample_dist(gen);

                sample_results = process_strats(
                    input_data,
                    params[p](p_pr_hosp_offset, t), params[p](p_los_offset, t),
                    i_prior, gen
                );

                arma::vec sampled_ward_occupancy(n_fit_days, arma::fill::zeros);
                arma::vec sampled_ICU_occupancy(n_fit_days, arma::fill::zeros);

                for(int s = 0; s < def_n_strats; s++) {
                    for(int x = day_start_fit * def_n_compartment_groups; x < day_start_fit + n_fit_days * def_n_compartment_groups; x++) {
                        int d = x % input_data.n_days;

                        if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ward) {
                            sampled_ward_occupancy(d - day_start_fit) += sample_results[s].grouped_occupancy_counts[x];
                        } else if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ICU) {
                            sampled_ICU_occupancy(d - day_start_fit) += sample_results[s].grouped_occupancy_counts[x];
                        }
                    }
                }

                float nrmse_ward = vec_normalised_rmse(known_ward_occupancy, sampled_ward_occupancy);
                float nrmse_ICU = vec_normalised_rmse(known_ICU_occupancy, sampled_ICU_occupancy);

                float weighted_error = nrmse_ward + 0.5 * nrmse_ICU;

                if(weighted_error < input_data.thresholds[t])
                    fit_achieved = true;
				
				n_attempts++;
            }
            
			attempts(p, t) = n_attempts;

            if(t == input_data.thresholds.size() - 1) {
                results[p] = sample_results;

            } else if (t == 0) {
				weights(p, t) = 1;
			} else {
				double numer = 1;
				for(int w = 0; w < n_params; w++)
					numer *= normal_pdf(params[p](w, t), 0, 1);
				
				double denom = 0;

				for(int j = 0; j < n_particles; j++) {
					double prob_product = 1;
					for(int w = 0; w < n_params; w++) 
						prob_product *= normal_pdf(params[j](w, t), params[j](w, t - 1), perturb_sigma);

					denom += weights(j, t - 1) * prob_product;
				}
				
				weights(p, t) = numer / denom;
			}
        });
    }

    smc_results results_obj;
    results_obj.results = results;
    results_obj.weights = weights;
    results_obj.attempts = attempts;
    

    return results_obj;
}



std::vector<mush_results> abc_smc::process_strats(
    const smc_inputs &input_data,

    float pr_hosp_offset,
    float los_offset,

    int i_prior,
    std::mt19937 &rng
) {
    std::vector<mush_results> sample_results(def_n_strats);

    mush_params sim_params {
        .n_days = input_data.n_days,
        .steps_per_day = input_data.steps_per_day,
        .n_delay_samples = input_data.n_delay_samples
    };

    for (int s = 0; s < def_n_strats; s++)
    {
        strat_data s_data = input_data.strat_datas[s][i_prior % input_data.n_strat_samples];


        std::vector<int> hospitalised_cases = sample_hospitalised_cases_const(
            input_data.case_curves[s][i_prior],
            input_data.pr_hosp_curves[s][i_prior],

            pr_hosp_offset,

            rng
        );

        sample_results[s] = musher::mush_curve(
            sim_params,
            hospitalised_cases,
            s_data,
            los_offset,
            input_data.pr_ICU_curves[s][i_prior]
        );
    }

    return sample_results;
}