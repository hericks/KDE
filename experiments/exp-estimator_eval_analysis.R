# Algorithms
loop_over_samples <- function(x, samples, f) {
  ret <- numeric(length(x))
  for (x0 in samples) {
    ret <- ret + f(x - x0)
  }
  ret/length(samples)
}

loop_over_input <- function(x, samples, f) {
  ret <- numeric(length(x))
  for (i in seq_along(x)) {
    ret[i] <- sum(f(x[i] - samples))
  }
  ret/length(samples)
}

vec_rowMeans <- function(x, samples, f) {
  rowMeans(outer(x, samples, function(x, y) f(x - y)))
}

# Testing Definitions
optimal_algorithm <- function(num_x, num_samples, exprs_to_test, times = 20L) {
  print(paste0("num_x: ", num_x, ", num_samples: ", num_samples))
  mb <- microbenchmark::microbenchmark(list = exprs_to_test,
                                       times = times,
                                       setup = {
                                         x <- rnorm(num_x)
                                         samples <- rnorm(num_samples)
                                       })

  means <- tapply(mb$time, mb$expr, mean)
  return(names(means)[which.min(means)])
}

# Testing
exprs_to_test <- alist(
  LOS=loop_over_samples(x, samples, uniform_kernel_eval),
  LOI=loop_over_input(x, samples, uniform_kernel_eval),
  VEC=vec_rowMeans(x, samples, uniform_kernel_eval)
)

uniform_kernel_eval <- function(x) {
  ret <- double(length(x))
  ret[abs(x) <= 0.1] <- 1
  ret
}

num_x_to_test <- 2**(0:10)
num_samples_to_test <- 2**(0:10)

optimal_algorithms <- outer(num_x_to_test,
                            num_samples_to_test,
                            Vectorize(function(num_x, num_samples) {
                              optimal_algorithm(num_x, num_samples, exprs_to_test, times = 20L)
                            }))

dimnames(optimal_algorithms) <- list("num_x"=num_x_to_test, "num_samples"=num_samples_to_test)
optimal_algorithms
