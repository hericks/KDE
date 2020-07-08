#' Construct a kernel density estimator
#'
#' @description
#' The kernelDensityEstimator function factory takes observations, returning
#' the corresponding kernel density estimator.
#'
#' @param kernel The (vectorised) kernel function to use for the construction of the estimator satisfying [is_kernel].
#' @param samples A numerical vector to base the construction of the estimator on.
#' @param bandwidth A non-negative numeric value to use as the bandwidth for the kernel.
#'
#' @return The constructed kernel density estimator. Therefore a function vectorised in its evaluation point returning the estimator at these points.
#'
#'
#' @seealso
#' * [kernels()] for a list of already implemented kernels.
#' * [validate_kernel()] to validate cusotm kernel functions.
#'
#' @export
kernelDensityEstimator <- function(kernel, samples, bandwidth=1) {
  # Kernel conditions
  stopifnot("the kernel has to be valid" = validate_Kernel(kernel))

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
      ret <- ret + kernel_transform(kernel, x0, bandwidth)(x)
    }
    ret/length(samples)
  }
}

