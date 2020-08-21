#' Construct a Kernel Density Estimator
#'
#' @description The \code{kernel_density_estimator} function builds a kernel
#'   density estimator using the provided kernel, bandwidth and samples.
#'
#' @param kernel a kernel as S3 object of the class \link{Kernel}.
#' @param samples numeric vector; the observations.
#' @param bandwidth non-negative numeric scalar; the bandwidth for the
#'   estimator.
#' @param subdivisions positive numeric scalar; subdivisions parameter
#'   internally passed to \code{\link{integrate_primitive}}.
#'
#' @details The validation of the returned estimator as
#'   \code{\link{IntegrableFunction}} relies on the function
#'   \code{integrate_primitive}, thus the \code{subdivisions} parameter.
#'
#' @return The estimator as S3 object of class \code{\link{IntegrableFunction}}.
#'
#' @source Nonparametric Estimation, Comte \[2017\], ISBN: 978-2-36693-030-6
#'
#' @seealso \code{\link{Kernel}} for more information about kernels,
#'   \code{\link{pco_method}}, \code{\link{cross_validation}} and
#'   \code{\link{goldenshluger_lepski}} for automatic bandwidth-selection
#'   algorithms.
#'
#' @export
kernel_density_estimator <- function(kernel, samples, bandwidth=1, subdivisions=1000L) {
  # Kernel conditions
  validate_Kernel(kernel)

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
