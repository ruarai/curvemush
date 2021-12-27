
#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame mush(
  int n_samples,
  int n_delay_samples,
  
  int n_days,
  int steps_per_day,
  
  int t_forecast_start,
  
  NumericMatrix ensemble_curves
) {
  
  mush_params params;
  params.n_days = n_days;
  params.steps_per_day = steps_per_day;
  params.n_delay_samples = n_delay_samples;
  
  std::vector<std::vector<int>> hosp_curves(n_samples);
  
  int n_curves = ensemble_curves.ncol();
  
  for(int i = 0; i < n_samples; i++) {
    NumericVector curve_i = ensemble_curves(_, i % n_curves);
    
    hosp_curves[i] = std::vector<int>(n_days, -1);
    
    for(int d = 0; d < t_forecast_start; d++) {
      hosp_curves[i][d] = curve_i[d * 9];
    }
    
    for(int d = 0; d < n_days - t_forecast_start; d++) {
      int n_cases = curve_i[t_forecast_start * 9 + d];
      
      for(int j = 0; j < n_cases; j++) {
        if(R::runif(0, 1) < 0.1)
          hosp_curves[i][t_forecast_start + d]++;
      }

    }
  }
  
  std::vector<mush_results> results(n_samples);
  RcppThread::parallelFor(
    0, n_samples,
    [&results, params, &hosp_curves] (unsigned int i) {
    results[i] = musher::mush_curve(params, hosp_curves[i]);
  });
  
  int result_size = params.n_days * def_n_compartment_groups;
  int n_results = n_samples * result_size;
  
  NumericVector sample_label(n_results);
  NumericVector t_day(n_results);
  NumericVector compartment_group_label(n_results);
  NumericVector compartment_group_count(n_results);
  
  for(int i = 0; i < n_samples; i++) {
    for(int x = 0; x < result_size; x++) {
      int ix = i * result_size + x;
      
      sample_label[ix] = i;
      t_day[ix] = x % params.n_days;
      compartment_group_label[ix] = results[i].grouped_occupancy_compartment_labels[x];
      compartment_group_count[ix] = results[i].grouped_occupancy_counts[x];
    }
  }
  
  DataFrame R_results = DataFrame::create(
    _["sample"] = sample_label,
    _["t_day"] = t_day,
    _["compartment_group"] = compartment_group_label,
    _["count"] = compartment_group_count
  );
  
  
  return R_results;
}
