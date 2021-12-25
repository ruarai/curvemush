library(curvemush)

n_compartments <- 4

ix <- function(t, compartment, type) {
  return((t - 1) * n_compartments * 2 + (compartment - 1) * 2 + type)
}

n_many <- 10000
a <- Sys.time()

results <- mush(n_many,
                1000,
                512)

b <- Sys.time()


n_time <- length(results[[1]]) / (2 * n_compartments * 4)
n_time <- 40

plot(results[[1]][ix(1:n_time, 1, 1)], type ='l', col = rgb(1,0,0,0.1), ylim = c(0,n_many))

for(j in 2:100) {
  lines(results[[j]][ix(1:n_time, 1, 1)], type ='l', col= rgb(1,0,0,0.1))
  
  for(i in 2:n_compartments){
    p <- i/n_compartments
    lines(results[[j]][ix(1:n_time, i, 1)], type ='l', col = rgb(1-p,p/2,p,0.1))
  }
  
}

print(b - a)
