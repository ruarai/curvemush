
#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"
#include "mush_params.h"


// [[Rcpp::export]]
std::vector<std::vector<int>> mush(
  int n_samples,
  int n_delay_samples
) {
  
  std::vector<std::vector<int>> results(n_samples);
  
  mush_params params;
  params.n_compartments = 4;
  params.n_days = 365;
  params.steps_per_day = 4;
  params.n_delay_samples = n_delay_samples;
  
  
  RcppThread::parallelFor(0, n_samples,
                          [&results, params] (unsigned int i) {
    results[i] = musher::mush_curve(params);
  });
  
  return results;
}
