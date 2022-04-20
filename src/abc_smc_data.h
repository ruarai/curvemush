

#pragma once

#include <RcppArmadillo.h>
#include "mush_params.h"
#include "strat_data.h"
#include "musher.h"


class smc_inputs {
public:
    std::vector<int> known_ward_occupancy;
    std::vector<int> known_ICU_occupancy;

    std::vector<std::vector<std::vector<int>>> case_curves;
    std::vector<std::vector<std::vector<float>>> pr_hosp_curves;
    std::vector<std::vector<std::vector<float>>> pr_ICU_curves;

    std::vector<std::vector<strat_data>> strat_datas;

    std::vector<float> thresholds;

    arma::mat spline_basis;

    int n_days;
    int n_delay_samples;
    int steps_per_day;

    int n_particles;

    int n_parameter_samples;
    int n_strat_samples;

    int day_start_fit;
    int day_end_fit;
};

class smc_results {
public:
    std::vector<std::vector<mush_results>> results;
    std::vector<arma::vec> pr_hosp_fits;

    arma::mat weights;
    arma::mat attempts;
};