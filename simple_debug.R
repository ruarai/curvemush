


library(curvemush)

load(".debug")


a <- Sys.time()
results <- mush(
  n_samples = n_samples,
  n_delay_samples = n_delay_samples,
  
  n_days = n_days,
  steps_per_day = steps_per_day,
  
  t_forecast_start = as.numeric(forecast_start_date - sim_start),
  
  ensemble_curves = input_curves
)
b <- Sys.time()

print(b - a)