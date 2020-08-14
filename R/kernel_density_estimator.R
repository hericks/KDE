#' Construct a kernel density estimator
#'
#' @description The \code{kernel_density_estimator} function builds a kernel
#'   density estimator using the provided kernel, bandwidth and samples.
#'
#' @param kernel a kernel as S3 object of the class \link{Kernel}.
#' @param samples a numerical vector of observations.
#' @param bandwidth a non-negative numeric value to use as the bandwidth for the
#'   estimator.
#' @param subdivisions a integer vector of length 1 used for the subdivisions
#'   parameter of the builtin R-function \code{\link{integrate}}.
#'
#' @details The validation of the returned estimator as
#'   \code{\link{IntegrableFunction}} relies on the builtin function
#'   \code{\link{integrate}}, which requires a \code{subdivisions} argument.
#'   Integration using a larger number of subdivisions will increase the
#'   runtime. In contrast too few subdivisions may result in a runtime error.
#'
#'   For more information about kernel density estimaton, see "Nonparametric
#'   Estimation" by Fabienne Comte.
#'
#' @return The estimator is returned as S3 object of class
#'   \code{\link{IntegrableFunction}}.
#'
#' @source Nonparametric Estimation, Comte \[2017\], ISBN: 978-2-36693-030-6
#'
#' @seealso \code{\link{Kernel}} for more information about kernels,
#'   \code{\link[KDE:pco_method]{PCO}},
#'   \code{\link[KDE:cross_validation]{Cross-Validation}} and
#'   \code{\link[KDE:goldenshluger_lepski]{Goldenshluger-Lepski}} for
#'   automatic bandwidth-selection algorithms.
#'
#' @export
kernel_density_estimator <- function(kernel, samples, bandwidth=1, subdivisions=100L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  # Subdivisions conditions
  subdivisions <- subdivisions
  stopifnot("Entry 'subdivisions' must be numeric"=is.numeric(subdivisions))
  stopifnot("Entry 'subdivisions' must be positive"=subdivisions > 0)

  kernel_eval <- kernel$fun

  estimator_eval <- function(x) {
    stopifnot("input has to be numeric"=is.numeric(x))

    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel_eval((x - x0)/bandwidth)/bandwidth
    }
    ret/length(samples)
  }

  support <- bandwidth*kernel$support + range(samples)

  IntegrableFunction(estimator_eval, support, subdivisions=subdivisions)
}
