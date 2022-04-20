


library(curvemush)

load(".debug")


a <- Sys.time()
results <- curvemush::mush_abc(
  n_samples = 4000,
  n_delay_samples = 512,
  
  n_outputs = 10000,
  
  n_days = case_trajectories$n_days,
  steps_per_day = 4,
  
  thresholds_vec = thresholds,
  rejections_per_selections = 100,
  do_ABC = TRUE,
  
  prior_sigma_los = 0,
  prior_sigma_hosp = 0,
  
  ensemble_curves = case_curves,
  
  forecasting_parameters = forecasting_parameters,
  
  known_ward_vec = occupancy_curve_match$ward_vec,
  known_ICU_vec = occupancy_curve_match$ICU_vec,
  
  mat_pr_age_given_case = mat_pr_age_given_case,
  mat_pr_hosp = mat_pr_hosp,
  mat_pr_ICU = mat_pr_ICU
)

b <- Sys.time()


print(b-a)