#include "group_data.h"



group_data group_data::read_group_data(DataFrame forecasting_parameters, int group_ix) {

    group_data data;

    data.pr_age_given_case = as<NumericVector>(forecasting_parameters["pr_age_given_case"])[group_ix];
    data.pr_hosp = as<NumericVector>(forecasting_parameters["pr_hosp"])[group_ix];

    data.pr_ward_to_discharge = as<NumericVector>(forecasting_parameters["pr_ward_to_discharge"])[group_ix];
    data.pr_ward_to_ICU = as<NumericVector>(forecasting_parameters["pr_ward_to_ICU"])[group_ix];

    data.pr_ICU_to_discharge = as<NumericVector>(forecasting_parameters["pr_ICU_to_discharge"])[group_ix];
    data.pr_ICU_to_postICU = as<NumericVector>(forecasting_parameters["pr_ICU_to_postICU"])[group_ix];

    data.pr_postICU_to_death = as<NumericVector>(forecasting_parameters["pr_postICU_to_death"])[group_ix];

    data.d_shape_symptomatic_to_ward = as<NumericVector>(forecasting_parameters["shape_symptomatic_to_ED"])[group_ix];
    data.d_scale_symptomatic_to_ward = as<NumericVector>(forecasting_parameters["scale_symptomatic_to_ED"])[group_ix];
    
    data.d_shape_ward_to_discharge = as<NumericVector>(forecasting_parameters["shape_ward_to_discharge"])[group_ix];
    data.d_scale_ward_to_discharge = as<NumericVector>(forecasting_parameters["scale_ward_to_discharge"])[group_ix];
    
    data.d_shape_ward_to_ICU = as<NumericVector>(forecasting_parameters["shape_ward_to_ICU"])[group_ix];
    data.d_scale_ward_to_ICU = as<NumericVector>(forecasting_parameters["scale_ward_to_ICU"])[group_ix];
    
    data.d_shape_ward_to_death = as<NumericVector>(forecasting_parameters["shape_ward_to_death"])[group_ix];
    data.d_scale_ward_to_death = as<NumericVector>(forecasting_parameters["scale_ward_to_death"])[group_ix];
    
    data.d_shape_ICU_to_discharge = as<NumericVector>(forecasting_parameters["shape_ICU_to_discharge"])[group_ix];
    data.d_scale_ICU_to_discharge = as<NumericVector>(forecasting_parameters["scale_ICU_to_discharge"])[group_ix];
    
    data.d_shape_ICU_to_postICU = as<NumericVector>(forecasting_parameters["shape_ICU_to_postICU"])[group_ix];
    data.d_scale_ICU_to_postICU = as<NumericVector>(forecasting_parameters["scale_ICU_to_postICU"])[group_ix];

    data.d_shape_ICU_to_death = as<NumericVector>(forecasting_parameters["shape_ICU_to_death"])[group_ix];
    data.d_scale_ICU_to_death = as<NumericVector>(forecasting_parameters["scale_ICU_to_death"])[group_ix];
    
    data.d_shape_postICU_to_discharge = as<NumericVector>(forecasting_parameters["shape_postICU_to_discharge"])[group_ix];
    data.d_scale_postICU_to_discharge = as<NumericVector>(forecasting_parameters["scale_postICU_to_discharge"])[group_ix];

    data.d_shape_postICU_to_death = as<NumericVector>(forecasting_parameters["shape_postICU_to_death"])[group_ix];
    data.d_scale_postICU_to_death = as<NumericVector>(forecasting_parameters["scale_postICU_to_death"])[group_ix];

    return data;
}