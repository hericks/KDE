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

stack <- function(x, samples, f) {
  sorted_x <- sort(x)
  sorted_samples <- sort(samples)
  n1 <- integer(length(x))
  i1 <- 1
  i2 <- 1
  l <- length(x)
  for (x2 in sorted_samples) {
    # advance i1 until v1 is in range below
    while (x2-sorted_x[i1] > 0.001 & i1 <= l) i1 <- i1+1
    if (i2 > i1) n1[i1:(i2-1)] <- n1[i1:(i2-1)]+1 else i2 <- i1
    # advance i2 until out of range adding 1 to n1[i2] each time
    while (sorted_x[i2]-x2 <= 0.001 & i2 <= l) {
      n1[i2] <- n1[i2]+1
      i2 <- i2+1
    }
  }
  s5 <- n1[rank(x)]/length(samples)
  return(s5)
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
uniform_kernel_eval <- function(x) {
  ret <- double(length(x))
  ret[abs(x) <= 0.001] <- 1
  ret
}

exprs_to_test <- alist(
  LOS=loop_over_samples(x, samples, uniform_kernel_eval),
  LOI=loop_over_input(x, samples, uniform_kernel_eval),
  # VEC=vec_rowMeans(x, samples, uniform_kernel_eval),
  STK=stack(x, samples, uniform_kernel_eval)
)

num_x_to_test <- 2**(0:8)
num_samples_to_test <- 2**(0:8)

optimal_algorithms <- outer(num_x_to_test,
                            num_samples_to_test,
                            Vectorize(function(num_x, num_samples) {
                              optimal_algorithm(num_x, num_samples, exprs_to_test, times = 20L)
                            }))

dimnames(optimal_algorithms) <- list("num_x"=num_x_to_test, "num_samples"=num_samples_to_test)
optimal_algorithms


x <- rnorm(10000)
samples <- rnorm(1000)
microbenchmark::microbenchmark(loop_over_samples(x, samples, uniform_kernel_eval),
                               loop_over_input(x, samples, uniform_kernel_eval),
                               stack(x, samples, uniform_kernel_eval))
