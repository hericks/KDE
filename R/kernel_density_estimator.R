kernelDensityEstimator <- function(kernel, samples, bandwidth=1) {
  # Kernel conditions
  stopifnot("kernel has to satisfy is_kernel() conditions"=is_kernel(kernel))

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  force(kernel)
  force(samples)
  force(bandwidth)

  function(x) {
    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernelTransform(kernel, x0, bandwidth)(x)
    }
    ret/length(samples)
  }
}

