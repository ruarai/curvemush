// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mush_abc
List mush_abc(int n_samples, int n_delay_samples, int n_outputs, int n_days, int steps_per_day, NumericVector thresholds_vec, int rejections_per_selections, bool do_ABC, NumericMatrix ensemble_curves, NumericMatrix mat_pr_age_given_case, NumericMatrix mat_pr_hosp, NumericMatrix mat_pr_ICU, DataFrame forecasting_parameters, IntegerVector known_ward_vec, IntegerVector known_ICU_vec, double prior_sigma_los, double prior_sigma_hosp);
RcppExport SEXP _curvemush_mush_abc(SEXP n_samplesSEXP, SEXP n_delay_samplesSEXP, SEXP n_outputsSEXP, SEXP n_daysSEXP, SEXP steps_per_daySEXP, SEXP thresholds_vecSEXP, SEXP rejections_per_selectionsSEXP, SEXP do_ABCSEXP, SEXP ensemble_curvesSEXP, SEXP mat_pr_age_given_caseSEXP, SEXP mat_pr_hospSEXP, SEXP mat_pr_ICUSEXP, SEXP forecasting_parametersSEXP, SEXP known_ward_vecSEXP, SEXP known_ICU_vecSEXP, SEXP prior_sigma_losSEXP, SEXP prior_sigma_hospSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_delay_samples(n_delay_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_outputs(n_outputsSEXP);
    Rcpp::traits::input_parameter< int >::type n_days(n_daysSEXP);
    Rcpp::traits::input_parameter< int >::type steps_per_day(steps_per_daySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds_vec(thresholds_vecSEXP);
    Rcpp::traits::input_parameter< int >::type rejections_per_selections(rejections_per_selectionsSEXP);
    Rcpp::traits::input_parameter< bool >::type do_ABC(do_ABCSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ensemble_curves(ensemble_curvesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_age_given_case(mat_pr_age_given_caseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_hosp(mat_pr_hospSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_ICU(mat_pr_ICUSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type forecasting_parameters(forecasting_parametersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type known_ward_vec(known_ward_vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type known_ICU_vec(known_ICU_vecSEXP);
    Rcpp::traits::input_parameter< double >::type prior_sigma_los(prior_sigma_losSEXP);
    Rcpp::traits::input_parameter< double >::type prior_sigma_hosp(prior_sigma_hospSEXP);
    rcpp_result_gen = Rcpp::wrap(mush_abc(n_samples, n_delay_samples, n_outputs, n_days, steps_per_day, thresholds_vec, rejections_per_selections, do_ABC, ensemble_curves, mat_pr_age_given_case, mat_pr_hosp, mat_pr_ICU, forecasting_parameters, known_ward_vec, known_ICU_vec, prior_sigma_los, prior_sigma_hosp));
    return rcpp_result_gen;
END_RCPP
}
// mush_abc_smc
List mush_abc_smc(int n_parameter_samples, int n_particles, int n_delay_samples, int n_days, int steps_per_day, NumericMatrix ensemble_curves, NumericMatrix mat_pr_age_given_case, NumericMatrix mat_pr_hosp, NumericMatrix mat_pr_ICU, DataFrame forecasting_parameters, std::vector<float> thresholds_vec, std::vector<int> known_ward_vec, std::vector<int> known_ICU_vec);
RcppExport SEXP _curvemush_mush_abc_smc(SEXP n_parameter_samplesSEXP, SEXP n_particlesSEXP, SEXP n_delay_samplesSEXP, SEXP n_daysSEXP, SEXP steps_per_daySEXP, SEXP ensemble_curvesSEXP, SEXP mat_pr_age_given_caseSEXP, SEXP mat_pr_hospSEXP, SEXP mat_pr_ICUSEXP, SEXP forecasting_parametersSEXP, SEXP thresholds_vecSEXP, SEXP known_ward_vecSEXP, SEXP known_ICU_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_parameter_samples(n_parameter_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_particles(n_particlesSEXP);
    Rcpp::traits::input_parameter< int >::type n_delay_samples(n_delay_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_days(n_daysSEXP);
    Rcpp::traits::input_parameter< int >::type steps_per_day(steps_per_daySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ensemble_curves(ensemble_curvesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_age_given_case(mat_pr_age_given_caseSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_hosp(mat_pr_hospSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_pr_ICU(mat_pr_ICUSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type forecasting_parameters(forecasting_parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type thresholds_vec(thresholds_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type known_ward_vec(known_ward_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type known_ICU_vec(known_ICU_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(mush_abc_smc(n_parameter_samples, n_particles, n_delay_samples, n_days, steps_per_day, ensemble_curves, mat_pr_age_given_case, mat_pr_hosp, mat_pr_ICU, forecasting_parameters, thresholds_vec, known_ward_vec, known_ICU_vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_curvemush_mush_abc", (DL_FUNC) &_curvemush_mush_abc, 17},
    {"_curvemush_mush_abc_smc", (DL_FUNC) &_curvemush_mush_abc_smc, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_curvemush(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
