<<<<<<< HEAD
library(rlang)

kde <- function(y, h, kernel){
  n <- length(y)
  ker <- eval(expr(!!kernel))
  function(x){
    sapply(x,function(x) (1/(n*h))*sum(sapply(y, function(y) ker((x-y)/h))))
=======
kernelDensityEstimator <- function(kernel, samples, bandwidth=1) {
  # TODO: Condition checking for arguments

  force(kernel)
  force(samples)

  function(x) {
    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel((x - x0)/bandwidth)/bandwidth
    }
    ret/length(samples)
>>>>>>> 7132fe7005b9929b132653fd8db31e6693185551
  }
}
