#include "strat_data.h"



strat_data strat_data::read_strat_data(DataFrame forecasting_parameters, int strat_ix) {

    strat_data data;

    data.pr_ward_to_discharge = as<NumericVector>(forecasting_parameters["pr_ward_to_discharge"])[strat_ix];

    data.pr_ICU_to_discharge = as<NumericVector>(forecasting_parameters["pr_ICU_to_discharge"])[strat_ix];
    data.pr_ICU_to_postICU = as<NumericVector>(forecasting_parameters["pr_ICU_to_postICU"])[strat_ix];

    data.pr_postICU_to_death = as<NumericVector>(forecasting_parameters["pr_postICU_to_death"])[strat_ix];

    data.d_shape_symptomatic_to_ward = as<NumericVector>(forecasting_parameters["shape_onset_to_ward"])[strat_ix];
    data.d_scale_symptomatic_to_ward = as<NumericVector>(forecasting_parameters["scale_onset_to_ward"])[strat_ix];
    
    data.d_shape_ward_to_discharge = as<NumericVector>(forecasting_parameters["shape_ward_to_discharge"])[strat_ix];
    data.d_scale_ward_to_discharge = as<NumericVector>(forecasting_parameters["scale_ward_to_discharge"])[strat_ix];
    
    data.d_shape_ward_to_ICU = as<NumericVector>(forecasting_parameters["shape_ward_to_ICU"])[strat_ix];
    data.d_scale_ward_to_ICU = as<NumericVector>(forecasting_parameters["scale_ward_to_ICU"])[strat_ix];
    
    data.d_shape_ward_to_death = as<NumericVector>(forecasting_parameters["shape_ward_to_death"])[strat_ix];
    data.d_scale_ward_to_death = as<NumericVector>(forecasting_parameters["scale_ward_to_death"])[strat_ix];
    
    data.d_shape_ICU_to_discharge = as<NumericVector>(forecasting_parameters["shape_ICU_to_discharge"])[strat_ix];
    data.d_scale_ICU_to_discharge = as<NumericVector>(forecasting_parameters["scale_ICU_to_discharge"])[strat_ix];
    
    data.d_shape_ICU_to_postICU = as<NumericVector>(forecasting_parameters["shape_ICU_to_postICU"])[strat_ix];
    data.d_scale_ICU_to_postICU = as<NumericVector>(forecasting_parameters["scale_ICU_to_postICU"])[strat_ix];

    data.d_shape_ICU_to_death = as<NumericVector>(forecasting_parameters["shape_ICU_to_death"])[strat_ix];
    data.d_scale_ICU_to_death = as<NumericVector>(forecasting_parameters["scale_ICU_to_death"])[strat_ix];
    
    data.d_shape_postICU_to_discharge = as<NumericVector>(forecasting_parameters["shape_postICU_to_discharge"])[strat_ix];
    data.d_scale_postICU_to_discharge = as<NumericVector>(forecasting_parameters["scale_postICU_to_discharge"])[strat_ix];

    data.d_shape_postICU_to_death = as<NumericVector>(forecasting_parameters["shape_postICU_to_death"])[strat_ix];
    data.d_scale_postICU_to_death = as<NumericVector>(forecasting_parameters["scale_postICU_to_death"])[strat_ix];

    return data;
}