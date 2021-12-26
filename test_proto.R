library(curvemush)

library(ggplot2)





n_samples <- 10000
n_many <- 1000
steps_per_day <- 8
a <- Sys.time()

results <- mush(n_samples,
                512,
                steps_per_day)

b <- Sys.time()

compartment_labels <- c(
  "susceptible", "symptomatic", "ward", "discharged_ward", "died_ward", "ICU", "discharged_ICU", "died_ICU",
  "postICU_to_discharge", "postICU_to_death", "discharged_postICU", "died_postICU"
)

results$compartment <- compartment_labels[results$compartment + 1]



ggplot(results %>% filter(t_step == max(t_step), sample < 1000)) +
  geom_point(aes(y = compartment, x = count))


ggplot(results %>% filter(sample < 100)) +
  geom_line(aes(x = t_step, y = count, group = interaction(compartment, sample), color = compartment))


print(b - a)
