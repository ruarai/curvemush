library(curvemush)

n_compartments <- 5

ix <- function(t, compartment, type) {
  return((t - 1) * n_compartments * 2 + (compartment - 1) * 2 + type)
}

n_samples <- 1000
n_many <- 1000
steps_per_day <- 4
a <- Sys.time()

results <- mush(n_samples,
                512,
                steps_per_day)

b <- Sys.time()


n_time <- length(results[[1]]) / (2 * n_compartments * steps_per_day)

plot(results[[1]][ix(1:n_time, 1, 1)], type ='l', col = rgb(1,0,0,0.1), ylim = c(0,n_many))

for(j in 1:100) {
  lines(results[[j]][ix(1:n_time, 1, 1)], type ='l', col= rgb(1,0,0,0.1))
  
  for(i in 2:n_compartments){
    p <- i/n_compartments
    lines(results[[j]][ix(1:n_time, i, 1)], type ='l', col = rgb(1-p,p/2,p,0.1))
  }
  
}

print(b - a)
