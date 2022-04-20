


#pragma once

#include <RcppArmadillo.h>
using namespace Rcpp;

IntegerVector rmultinom_vec(unsigned int size, NumericVector &probs);