
#include <random>

#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"
#include "rmultinom_vec.h"

using namespace Rcpp;


// [[Rcpp::export]]
List mush(
  int n_samples,
  int n_delay_samples,
  
  int n_days,
  int steps_per_day,
  
  int t_forecast_start,
  
  NumericMatrix ensemble_curves,

  DataFrame forecasting_parameters,


  float scale_los
) {
  mush_params params;
  params.n_days = n_days;
  params.steps_per_day = steps_per_day;
  params.n_delay_samples = n_delay_samples;


  int n_strat_samples = forecasting_parameters.nrow() / def_n_strats;
  std::vector<strat_data> strat_datas(n_samples * def_n_strats);

  for(int i = 0; i < n_samples * def_n_strats; i++)
    strat_datas[i] = strat_data::read_strat_data(forecasting_parameters, i);
  

  // Not currently resampled
  NumericVector pr_age_given_case(def_n_strats);
  NumericVector pr_hosp_given_case(def_n_strats);
  for(int s = 0; s < def_n_strats; s++) 
  {
    pr_age_given_case[s] = strat_datas[s].pr_age_given_case;
    pr_hosp_given_case[s] = strat_datas[s].pr_hosp;
  }
  

  // Define and initialize our hospitalisation curves (vector of vectors of vectors)
  std::vector<std::vector<std::vector<int>>> hosp_curves(
    def_n_strats,
    std::vector<std::vector<int>>(n_samples, std::vector<int>(n_days, 0))
  );

  int n_curves = ensemble_curves.ncol();
  // Producing n_samples case/hospitalisation trajectories to feed into the model:
  for(int i = 0; i < n_samples; i++) {
    // Select a curve from our backcast/ensemble curve matrix
    // Looping (by % n_curves) if we are performing more samples than we have ensembles
    NumericVector curve_i = ensemble_curves(_, i % n_curves);
    
    
    // Over the backcast period, cases are split by age group (and known precisely)
    for(int s = 0; s < def_n_strats; s++)
      for(int d = 0; d < t_forecast_start; d++)
        hosp_curves[s][i][d] = curve_i[d * def_n_strats + s];
    
    // But in the forecast period, we perform the age sampling ourselves
    for(int d = 0; d < n_days - t_forecast_start; d++) {
      int n_cases = curve_i[t_forecast_start * def_n_strats + d];

      if(n_cases == 0)
        continue;

      // Produce a vector distributing the n_cases across each age stratification 
      // according to pr_age_given_case
      IntegerVector age_samples = rmultinom_vec(n_cases, pr_age_given_case);

      // Take our sampled cases and sample again according to pr_hosp for the age stratification
      for(int s = 0; s < def_n_strats; s++) {
        int n_hospitalised = R::rbinom(age_samples[s], pr_hosp_given_case[s]);

        hosp_curves[s][i][t_forecast_start + d] = n_hospitalised;
      }
    }
  }

  std::vector<std::vector<mush_results>> results(def_n_strats);
  
  for(int s = 0; s < def_n_strats; s++) {


    results[s] = std::vector<mush_results>(n_samples);
  
    // Perform our simulations in parallel using RcppThread
    // Important: can't use R data structures within threads!
    RcppThread::parallelFor(
      0, n_samples,
      [&results, params, &hosp_curves, s, strat_datas, n_strat_samples, n_days,
       scale_los] (unsigned int i) {

        strat_data s_data = strat_datas[i % n_strat_samples + s];


        results[s][i] = musher::mush_curve(params, hosp_curves[s][i], s_data, scale_los, std::vector<float>(n_days, s_data.pr_ward_to_ICU));
    }); 

  }


  int result_size = params.n_days * def_n_compartments;
  int n_results = n_samples * result_size;

  NumericVector sample_label(n_results);
  NumericVector t_day(n_results);
  NumericVector compartment_label(n_results);
  NumericVector compartment_count(n_results);
  NumericVector compartment_transitions(n_results);

  
  for(int i = 0; i < n_samples; i++) {
    for(int x = 0; x < result_size; x++) {
      int ix = i * result_size + x;
      
      sample_label[ix] = i;
      t_day[ix] = x % params.n_days;
      compartment_label[ix] = results[0][i].occupancy_compartment_labels[x];

      compartment_count[ix] = 0;
      compartment_transitions[ix] = 0;
      for(int s = 0; s < def_n_strats; s++) {
        compartment_count[ix] += results[s][i].occupancy_counts[x];
        compartment_transitions[ix] += results[s][i].transitions[x];
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
  int grped_n_results = n_samples * grped_result_size;
  
  NumericVector grped_sample_label(grped_n_results);
  NumericVector grped_t_day(grped_n_results);
  NumericVector grped_compartment_group_label(grped_n_results);
  NumericVector grped_compartment_group_count(grped_n_results);
  NumericVector grped_compartment_group_transitions(grped_n_results);
  
  for(int i = 0; i < n_samples; i++) {
    for(int x = 0; x < grped_result_size; x++) {
      int ix = i * grped_result_size + x;
      
      grped_sample_label[ix] = i;
      grped_t_day[ix] = x % params.n_days;
      grped_compartment_group_label[ix] = results[0][i].grouped_occupancy_compartment_labels[x];

      grped_compartment_group_count[ix] = 0;
      grped_compartment_group_transitions[ix] = 0;
      for(int s = 0; s < def_n_strats; s++) {
        grped_compartment_group_count[ix] += results[s][i].grouped_occupancy_counts[x];
        grped_compartment_group_transitions[ix] += results[s][i].grouped_transitions[x];
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


  List R_list_results = List::create(
    _["results"] = R_results,
    _["grouped_results"] = R_grped_results
  );
  
  
  return R_list_results;
}
