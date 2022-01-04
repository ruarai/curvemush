

setwd("/home/forecast/source/clinical_forecasting")

library(targets)

source("R/progression_model.R")


run_progression_model(
    tar_read(case_trajectories_NSW),
    tar_read(clinical_table_NSW),
    tar_read(nindss_state_NSW),

    tar_read(forecast_dates)
)