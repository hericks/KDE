kernelDensityEstimator <- function(kernel, samples, bandwidth=1) {
  # TODO: Condition checking for arguments
  force(kernel)
  force(samples)
  force(bandwidth)

  # kernel check
  stopifnot("kernel has to satisfy is_kernel() conditions"=is_kernel(kernel))
  # samples check
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)
  # bandwidth check
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  function(x) {
    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernelTransform(kernel, x0, bandwidth)
    }
    ret/length(samples)
  }
}

