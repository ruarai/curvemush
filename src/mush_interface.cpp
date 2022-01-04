
#include <random>

#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"

using namespace Rcpp;

IntegerVector rmultinom_vec(unsigned int size, NumericVector &probs) {
    int N = probs.size();

    IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}

// [[Rcpp::export]]
DataFrame mush(
  int n_samples,
  int n_delay_samples,
  
  int n_days,
  int steps_per_day,
  
  int t_forecast_start,
  
  NumericMatrix ensemble_curves,

  DataFrame forecasting_parameters
) {
  mush_params params;
  params.n_days = n_days;
  params.steps_per_day = steps_per_day;
  params.n_delay_samples = n_delay_samples;


  std::vector<group_data> g_datas(def_n_groups);

  for(int g = 0; g < def_n_groups; g++) 
    g_datas[g] = group_data::read_group_data(forecasting_parameters, g);

  NumericVector pr_age_given_case(def_n_groups);
  for(int g = 0; g < def_n_groups; g++) 
    pr_age_given_case[g] = g_datas[g].pr_age_given_case;
  

  // Define and initialize our hospitalisation curves
  std::vector<std::vector<std::vector<int>>> hosp_curves(def_n_groups);
  for(int g = 0; g < def_n_groups; g++) {
    hosp_curves[g] = std::vector<std::vector<int>>(n_samples);
    for(int i = 0; i < n_samples; i++)
      hosp_curves[g][i] = std::vector<int>(n_days, 0);
  }

  int n_curves = ensemble_curves.ncol();
  // Producing n_samples case/hospitalisation trajectories to feed into the model:
  for(int i = 0; i < n_samples; i++) {
    // Select a curve from our backcast/ensemble curve matrix
    // Looping (by % n_curves) if we are performing more samples than we have ensembles
    NumericVector curve_i = ensemble_curves(_, i % n_curves);
    
    
    // Over the backcast period, cases are split by age group (and known precisely)
    for(int g = 0; g < def_n_groups; g++)
      for(int d = 0; d < t_forecast_start; d++)
        hosp_curves[g][i][d] = curve_i[d * def_n_groups + g];
    
    // But in the forecast period, we perform the age sampling ourselves
    for(int d = 0; d < n_days - t_forecast_start; d++) {
      int n_cases = curve_i[t_forecast_start * def_n_groups + d];

      if(n_cases == 0)
        continue;

      IntegerVector age_samples = rmultinom_vec(n_cases, pr_age_given_case);

      for(int g = 0; g < def_n_groups; g++) {
        int n_cases_group = age_samples[g];

        int n_hospitalised = R::rbinom(n_cases_group, g_datas[g].pr_hosp);

        hosp_curves[g][i][t_forecast_start + d] = n_hospitalised;
      }
    }
  }

  std::vector<std::vector<mush_results>> results(def_n_groups);
  
  for(int g = 0; g < def_n_groups; g++) {

    group_data g_data = g_datas[g];

    results[g] = std::vector<mush_results>(n_samples);
  
    // Perform our simulations in parallel using RcppThread
    // Important: can't use R data structures within threads!
    RcppThread::parallelFor(
      0, n_samples,
      [&results, params, &hosp_curves, g, g_data] (unsigned int i) {
      results[g][i] = musher::mush_curve(params, hosp_curves[g][i], g_data);
    });

  }
  

  // Build out our result dataframe to return
  int result_size = params.n_days * def_n_compartment_groups;
  int n_results = n_samples * result_size;
  
  NumericVector sample_label(n_results);
  NumericVector t_day(n_results);
  NumericVector compartment_group_label(n_results);
  NumericVector compartment_group_count(n_results);
  NumericVector compartment_group_transitions(n_results);
  
  for(int i = 0; i < n_samples; i++) {
    for(int x = 0; x < result_size; x++) {
      int ix = i * result_size + x;
      
      sample_label[ix] = i;
      t_day[ix] = x % params.n_days;
      compartment_group_label[ix] = results[0][i].grouped_occupancy_compartment_labels[x];

      compartment_group_count[ix] = 0;
      compartment_group_transitions[ix] = 0;
      for(int g = 0; g < def_n_groups; g++) {
        compartment_group_count[ix] += results[g][i].grouped_occupancy_counts[x];
        compartment_group_transitions[ix] += results[g][i].grouped_transitions[x];
      }
    }
  }
  
  DataFrame R_results = DataFrame::create(
    _["sample"] = sample_label,
    _["t_day"] = t_day,
    _["compartment_group"] = compartment_group_label,
    _["count"] = compartment_group_count,
    _["transitions"] = compartment_group_transitions
  );
  
  
  return R_results;
}
