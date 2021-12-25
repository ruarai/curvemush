
library(curvemush)

n_samples <- 1000
n_many <- 100000
steps_per_day <- 1
a <- Sys.time()

results <- mush(n_samples,
                512,
                steps_per_day)

results[[1]]