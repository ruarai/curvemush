
#include <RcppArmadillo.h>
using namespace Rcpp;

// Helper functions to use rmultinom (multinomial sampling)
// Size is the sample size, probs is the probability density across [0, 1, ..., probs.size()]
IntegerVector rmultinom_vec(unsigned int size, NumericVector &probs) {
    int N = probs.size();

    IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}