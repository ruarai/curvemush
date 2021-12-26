
#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame mush(
  int n_samples,
  int n_delay_samples,
  int steps_per_day
) {
  
  std::vector<mush_results> results(n_samples);
  
  mush_params params;
  params.n_days = 100;
  params.steps_per_day = steps_per_day;
  params.n_delay_samples = n_delay_samples;
  
  
  RcppThread::parallelFor(0, n_samples,
                          [&results, params] (unsigned int i) {
    results[i] = musher::mush_curve(params);
  });
  
  int result_size = params.n_days * def_n_compartments;
  int n_results = n_samples * result_size;
  
  NumericVector sample_label(n_results);
  NumericVector t_day(n_results);
  NumericVector compartment_label(n_results);
  NumericVector compartment_count(n_results);
  
  for(int i = 0; i < n_samples; i++) {
    for(int x = 0; x < result_size; x++) {
      int ix = i * result_size + x;
      
      sample_label[ix] = i;
      t_day[ix] = x % params.n_days;
      compartment_label[ix] = results[i].occupancy_compartment_labels[x];
      compartment_count[ix] = results[i].occupancy_counts[x];
    }
  }
  
  DataFrame R_results = DataFrame::create(
    _["sample"] = sample_label,
    _["t_day"] = t_day,
    _["compartment"] = compartment_label,
    _["count"] = compartment_count
  );
  
  
  return R_results;
}
