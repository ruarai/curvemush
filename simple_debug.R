


library(curvemush)

load("../clinical_forecasting/.debug")


a <- Sys.time()
results <- curvemush::mush_abc(
  n_samples = 4000,
  n_delay_samples = 512,
  
  n_outputs = 1000,
  
  n_days = case_trajectories$n_days,
  steps_per_day = 4,
  
  thresholds_vec = thresholds,
  rejections_per_selections = 30,
  do_ABC = do_ABC,

  prior_sigma_los = prior_sigma_los,
  prior_sigma_hosp = prior_sigma_hosp,
  
  
  ensemble_curves = case_trajectories$curve_set,
  
  forecasting_parameters = clinical_parameter_samples,
  
  known_ward_vec = occupancy_curve_match$ward_vec,
  known_ICU_vec = occupancy_curve_match$ICU_vec,
  
  mat_pr_age_given_case = mat_pr_age_given_case,
  mat_pr_hosp = mat_pr_hosp,
  mat_pr_ICU = mat_pr_ICU
)

b <- Sys.time()


print(b-a)