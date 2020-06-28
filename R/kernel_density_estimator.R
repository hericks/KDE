kernelDensityEstimator <- function(kernel, samples, bandwidth=1) {
  # TODO: Condition checking for arguments

  force(kernel)
  force(samples)

  function(x) {
    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernelTransform(kernel, x0, bandwidth)
    }
    ret/length(samples)
  }
}
