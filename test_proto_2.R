
library(tidyverse)
library(lubridate)


assign_age_group <- function(age) {
  case_when(age < 10 ~ "0-9",
            age < 20 ~ "10-19",
            age < 30 ~ "20-29",
            age < 40 ~ "30-39",
            age < 50 ~ "40-49",
            age < 60 ~ "50-59",
            age < 70 ~ "60-69",
            age < 80 ~ "70-79",
            age >= 80 ~ "80+")
}

forecast_start_date <- ymd("2021-12-19")


nindss_curve <- read_rds("proto_data/nindss_nsw.rds") %>%
  filter(date_onset < forecast_start_date) %>%
  mutate(age_group = assign_age_group(age)) %>%
  #filter(ever_in_hospital) %>%
  
  group_by(date_onset, age_group) %>%
  summarise(n_hosp = n(), .groups = "drop") %>%
  complete(date_onset = full_seq(date_onset, 1), age_group) %>%
  
  mutate(n_hosp = replace_na(n_hosp, 0)) %>%
  arrange(date_onset, age_group)


# Need to match this to n columns once NA cols. removed below
nindss_curve_set <- nindss_curve$n_hosp %>%
  matrix(nrow = length(nindss_curve$n_hosp), ncol = 8000)


sim_start <- nindss_curve$date_onset %>% min()


ensemble_data <- read_csv("proto_data/combined_samples_delta_12021-12-18.csv")

ensemble_curves_df <- ensemble_data %>%
  filter(state == "NSW",
         date < min(ensemble_data$date) + ddays(28)) %>%
  select(-c(state, forecast_origin)) %>%

  pivot_wider(names_from = ".model",
              values_from = starts_with("sim")) %>%

  mutate(across(starts_with("sim"), ~ round(. * 0.1))) %>%

  arrange(date) %>%
  select(-date)


na_cols <- map_lgl(1:ncol(ensemble_curves_df), ~ all(is.na(ensemble_curves_df[,.])))

print(paste0("Dropping ", sum(na_cols), " columns for being entirely NA"))

ensemble_curves <- as.matrix(ensemble_curves_df[, !na_cols])

age_probs <- c(0.1, 0.1, 0.2, 0.2, 0, 0.1, 0.1, 0.15, 0.05)
input_curves <- rbind(nindss_curve_set, ensemble_curves)

n_backcast_days <- nrow(nindss_curve_set) / 9
n_days <- n_backcast_days + nrow(ensemble_curves)
steps_per_day <- 16

n_delay_samples <- 512

n_samples <- 2000



group_labels <- c("symptomatic", "ward", "ICU", "discharged", "died")

library(curvemush)

save.image(".debug")

a <- Sys.time()
results <- mush(
  n_samples = n_samples,
  n_delay_samples = n_delay_samples,
  
  n_days = n_days,
  steps_per_day = steps_per_day,
  
  t_forecast_start = n_backcast_days,
  
  ensemble_curves = input_curves
)
b <- Sys.time()

results_plot <- results %>%
  filter(sample < 100) %>%
  mutate(date = sim_start + ddays(t_day),
         group = group_labels[compartment_group + 1])


ggplot(results_plot) +
  geom_line(aes(x = date, y = count, group = sample)) +
  facet_wrap(~group, scales = "free_y") +
  
  coord_cartesian(c(forecast_start_date - 14, forecast_start_date + 28))


make_quants <- function(tbl) {
  data_matrix <- tbl %>%
    select(starts_with("sim_")) %>%
    as.matrix()
  
  id_tbl <- tbl %>%
    select(!starts_with("sim_"))
  
  medians <- data_matrix %>%
    matrixStats::rowMedians() %>%
    tibble(median = .)
  
  probs <- c(0.5, 0.75, 0.9, 0.95, 0.99)
  
  quant_probs <- c(rev(1 - probs) / 2, 0.5 + probs / 2)
  quant_names <- c(str_c("lower_", rev(probs) * 100), str_c("upper_", probs * 100))
  
  quants <- data_matrix %>%
    matrixStats::rowQuantiles(probs = quant_probs) %>%
    `colnames<-`(quant_names) %>%
    as_tibble() %>%
    bind_cols(id_tbl, .) %>%
    pivot_longer(cols = -all_of(colnames(id_tbl)),
                 names_to = c("type", "quant"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "type",
                values_from = "value") %>%
    
    mutate(quant = factor(quant, levels = as.character(probs * 100)) %>% fct_rev())
  
  quants
}

results_quants <- results %>%
  pivot_wider(names_from = "sample",
              names_prefix = "sim_",
              values_from = "count") %>%
  make_quants() %>%
  mutate(date = sim_start + ddays(t_day),
         group = group_labels[compartment_group + 1])


ggplot(results_quants %>% 
         filter(group == "ward" | group == "ICU")) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, group = quant),
              fill = 'blue3', alpha = 0.2) +
  
  facet_wrap(~group,
             scales = "free_y")

print(b - a)
