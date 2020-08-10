#' Construct a kernel density estimator
#'
#' @description The kernelDensityEstimator function factory takes observations,
#' returning the corresponding kernel density estimator.
#'
#' @param kernel The (vectorised) kernel function to use for the construction of
#'   the estimator satisfying [is_kernel].
#' @param samples A numerical vector to base the construction of the estimator
#'   on.
#' @param bandwidth A non-negative numeric value to use as the bandwidth for the
#'   kernel.
#'
#' @return The constructed kernel density estimator. Therefore a function
#'   vectorised in its evaluation point returning the estimator at these points.
#'
#' @seealso * [kernels()] for a list of already implemented kernels. *
#' [validate_kernel()] to validate cusotm kernel functions.
#'
#' @export
kernelDensityEstimator <- function(kernel, samples, bandwidth = 1, subdivisions = 100L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  kernel_eval <- kernel$fun

  estimator_eval <- function(x) {
    stopifnot("input has to be numeric"=is.numeric(x))

    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel_eval((x - x0)/bandwidth)
    }
    ret/(bandwidth*length(samples))
  }

  support <- c(bandwidth*kernel$support[1] + min(samples), bandwidth*kernel$support[2] + max(samples))

  IntegrableFunction(estimator_eval, support, subdivisions)
}







