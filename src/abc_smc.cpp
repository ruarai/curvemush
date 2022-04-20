#include "abc_smc.h"
#include "abc_smc_data.h"

#include "hosp_sampling.h"

#include <RcppArmadillo.h>
#include <RcppThread.h>
using namespace Rcpp;


typedef std::vector<double> stdvec;


double normal_pdf(double x, double m, double s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

smc_results abc_smc::process_abc_smc(
    smc_inputs &input_data
) {

    int n_particles = input_data.n_particles;
    int n_thresholds = input_data.thresholds.size();

    int df_spline = input_data.spline_basis.n_cols;

    std::vector<arma::mat> params(n_particles, arma::mat(df_spline, n_thresholds, arma::fill::zeros));

    arma::mat weights(n_particles, n_thresholds);
    arma::mat attempts(n_particles, n_thresholds);

    
    std::vector<std::vector<mush_results>> results(n_particles);

    std::vector<arma::vec> pr_hosp_fits(n_particles);

    double perturb_sigma = 0.05;

    for(int t = 0; t < input_data.thresholds.size(); t++) {
		Rcout << "t = " << t << ", (" << input_data.thresholds[t] << ")\n";

        RcppThread::parallelFor(0, n_particles,
        
        [&params, &weights, &attempts, &results, &pr_hosp_fits, &input_data,
         df_spline, perturb_sigma, n_particles, t] (unsigned int p) {
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
					for (int i = 0; i < df_spline; i++)
						params[p](i, t) = prior(gen);
				}
				else {
					std::vector<double> prev_weights = arma::conv_to<stdvec>::from(weights.col(t - 1));
					std::discrete_distribution<> sample_dist(prev_weights.begin(), prev_weights.end());

					int p_sample = sample_dist(gen);

					for(int w = 0; w < df_spline; w++)
						params[p](w, t) = params[p_sample](w, t - 1) + perturb_kernel(gen);
					
				}

                spline_params = params[p].col(t);
                
                
                i_prior = param_sample_dist(gen);

                sample_results = process_strats(
                    input_data, i_prior, spline_params, gen
                );
                std::vector<int> sample_ward_counts(input_data.n_days, 0);
                std::vector<int> sample_ICU_counts(input_data.n_days, 0);

                for(int s = 0; s < def_n_strats; s++) {
                    for(int x = 0; x < input_data.n_days * def_n_compartment_groups; x++) {

                        if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ward) {
                            sample_ward_counts[x % input_data.n_days] += sample_results[s].grouped_occupancy_counts[x];
                        } else if(sample_results[s].grouped_occupancy_compartment_labels[x] == cgrp_ICU) {
                            sample_ICU_counts[x % input_data.n_days] += sample_results[s].grouped_occupancy_counts[x];
                        }
                    }
                }

                
                bool rejected = false;

                for(int d = 0; d < input_data.n_days; d += 4) {
                    if(input_data.known_ward_occupancy[d] != -1 && 
                        std::abs(sample_ward_counts[d] - input_data.known_ward_occupancy[d]) > std::max(input_data.known_ward_occupancy[d] * input_data.thresholds[t], 10.0f)) {
                        rejected = true;
                    }
                    // if(known_ICU[t] != -1 &&
                    //     std::abs(sample_ICU_counts[t] - known_ICU[t]) > std::max(known_ICU[t] * thresholds[i_threshold] * 2, 10.0f)) {
                    //     rejected = true;
                    // }
                }

                if(!rejected)
                    fit_achieved = true;
				
				n_attempts++;
            }
            
			attempts(p, t) = n_attempts;

            if(t == input_data.thresholds.size() - 1) {
                results[p] = sample_results;

                pr_hosp_fits[p] = input_data.spline_basis * spline_params;

            } else if (t == 0) {
				weights(p, t) = 1;
			} else {
				double numer = 1;
				for(int w = 0; w < df_spline; w++)
					numer *= normal_pdf(params[p](w, t), 0, 1);
				
				double denom = 0;

				for(int j = 0; j < n_particles; j++) {
					double prob_product = 1;
					for(int w = 0; w < df_spline; w++) 
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
    results_obj.pr_hosp_fits = pr_hosp_fits;
    

    return results_obj;
}



std::vector<mush_results> abc_smc::process_strats(
    const smc_inputs &input_data,

    int i_prior,
    arma::vec spline_params,
    std::mt19937 &rng
) {
    std::vector<mush_results> sample_results(def_n_strats);

    arma::vec pr_hosp_offset = input_data.spline_basis * spline_params;


    mush_params sim_params {
        .n_days = input_data.n_days,
        .steps_per_day = input_data.steps_per_day,
        .n_delay_samples = input_data.n_delay_samples
    };

    for (int s = 0; s < def_n_strats; s++)
    {
        strat_data s_data = input_data.strat_datas[s][i_prior % input_data.n_strat_samples];


        std::vector<int> hospitalised_cases = sample_hospitalised_cases_spline_B(
            input_data.case_curves[s][i_prior],
            input_data.pr_hosp_curves[s][i_prior],

            pr_hosp_offset,
            input_data.day_start_fit, input_data.day_end_fit,

            rng
        );


        
        sample_results[s] = musher::mush_curve(
            sim_params,
            hospitalised_cases,
            s_data,
            0,
            input_data.pr_ICU_curves[s][i_prior]
        );
    }

    return sample_results;
}