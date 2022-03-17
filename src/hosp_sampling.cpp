
#include <Rcpp.h>
#include <random>

float logit(float x) { return std::log(x / (1 - x)); }
float ilogit(float x) { return std::exp(x) / (std::exp(x) + 1);}

std::vector<int> sample_hospitalised_cases(
    std::vector<int> &daily_cases,
    std::vector<float> &pr_hosp,
    float pr_hosp_scale,

    std::mt19937 &rng
) {
    int n_days = daily_cases.size();

    std::vector<int> daily_hospitalised(n_days, 0);

    for (int d = 0; d < n_days; d++)
    {
        int n_cases = daily_cases[d];

        if (n_cases == 0)
            continue;

        float adj_pr_hosp = ilogit(logit(pr_hosp[d]) + pr_hosp_scale);
            
        std::binomial_distribution<int> hosp_binom(n_cases, adj_pr_hosp);

        daily_hospitalised[d] = hosp_binom(rng);
    }

    return daily_hospitalised;
}