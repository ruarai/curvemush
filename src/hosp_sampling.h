

#pragma once
#include <RcppArmadillo.h>

#include <random>


std::vector<int> sample_hospitalised_cases(
    std::vector<int> &daily_cases,
    std::vector<float> &pr_hosp,
    float pr_hosp_scale,

    std::mt19937 &rng
);


std::vector<int> sample_hospitalised_cases_const(
    const std::vector<int> &daily_cases,
    const std::vector<float> &pr_hosp,
    float pr_hosp_scale,

    std::mt19937 &rng
);
