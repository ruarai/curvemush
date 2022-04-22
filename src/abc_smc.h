
#pragma once

#include "abc_smc_data.h"
#include "musher.h"
#include <RcppArmadillo.h>
using namespace Rcpp;


class abc_smc {
public:
    static smc_results process_abc_smc(
        smc_inputs &input_data
    );

    static std::vector<mush_results> process_strats(
        const smc_inputs &input_data,
        float pr_hosp_offset,
        float los_offset,
        int i_prior,
        std::mt19937 &rng
    );
};

