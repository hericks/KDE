x <- rnorm(250)
samples <- rnorm(250)
kernel <- rectangular
bandwidth <- 1

f1 <- function(x, samples, kernel, bandwidth) {
  ret <- numeric(length(x))
  for (x0 in samples) {
    ret <- ret + kernel$fun((x - x0)/bandwidth)
  }
  ret/(bandwidth*length(samples))
}

f2 <- function(x, samples, kernel, bandwidth) {
  apply(kernel$fun(outer(x, samples, `-`)/bandwidth), 1, sum)/(bandwidth*length(samples))
}

f1(x, samples, kernel, bandwidth)
f2(x, samples, kernel, bandwidth)

library(microbenchmark)
microbenchmark(f1(x, samples, kernel, bandwidth), f2(x, samples, kernel, bandwidth), times=10)
