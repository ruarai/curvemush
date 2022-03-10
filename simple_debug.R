


library(curvemush)

load(".debug")


a <- Sys.time()
results <- curvemush::mush_abc(
  n_samples = 400,
  n_delay_samples = 512,
  
  n_outputs = 100,
  
  n_days = case_trajectories$n_days,
  steps_per_day = 16,
  
  ward_threshold = 1000,
  
  prior_sigma_los = prior_sigma_los,
  prior_sigma_hosp = prior_sigma_hosp,
  
  t_forecast_start = case_trajectories$step_sampling_start,
  
  ensemble_curves = case_curves,
  
  forecasting_parameters = forecasting_parameters,
  
  known_ward_vec = occupancy_curve_match$count_vec,
  
  mat_pr_age_given_case = mat_pr_age_given_case,
  mat_pr_hosp = mat_pr_hosp,
  mat_pr_ICU = mat_pr_ICU
)
b <- Sys.time()


print(b-a)