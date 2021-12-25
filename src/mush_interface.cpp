
#include <Rcpp.h>
#include <RcppThread.h>
#include "musher.h"


// [[Rcpp::export]]
std::vector<std::vector<int>> mush(
  int n_start,
  
  int n_samples
) {
  
  std::vector<std::vector<int>> results(n_samples);
  
  
  RcppThread::parallelFor(0, n_samples, [&results, n_start] (unsigned int i) {
    results[i] = musher::mush_curve(n_start);
  });
  
  return results;
}
