#' Cross-Validation
#'
#' @description The cross-validation method is used to estimate an optimal
#'   bandwidth for kernel density estimation from a given set of bandwidths.
#'
#' @param kernel S3 object of class \code{\link{Kernel}}; the kernel to use for
#'   the estimator
#' @param samples numeric vector; the observations.
#' @param bandwidths strictly positive numeric vector; the bandwidth set from
#'   which the bandwidth with the least estimated risk will be selected.
#' @param subdivisions positive numeric scalar; subdivisions parameter
#'   internally passed to \code{\link{integrate_primitive}}.
#'
#' @details Cross-validation aims to minimize the mean integrated squared error
#'   (MISE) of a kernel density estimator. The MISE is defined as the
#'   expectation of the squared L2-Norm of the difference between estimator and
#'   (unknown) true density.
#'
#'   For each bandwidth \code{h} given in \code{bandwidths},
#'   \code{cross_validation} approximates the estimator-dependent part of the
#'   risk. The method then selects the bandwidth with the minimal associated
#'   risk.
#'
#' @return The estimated optimal bandwidth contained in the bandwidth set.
#'
#' @seealso \code{\link{kernel_density_estimator}} for more information about
#'   kernel density estimators, \code{\link{pco_method}} and
#'   \code{\link{goldenshluger_lepski}} for more automatic bandwidth-selection
#'   algorithms.
#'
#' @source \href{https://spartacus-idh.com/030.html}{Nonparametric Estimation},
#'   Comte \[2017\], ISBN: 978-2-36693-030-6
#'
#' @include kernel.R
#' @include kernel_density_estimator.R
#' @include logarithmic_bandwidth_set.R
#'
#' @export
cross_validation <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), subdivisions = 1000L) {
  # conditions for kernel
  tryCatch({
    validate_Kernel(kernel)
  }, error = "the kernel has to be valid")

  # conditions for samples
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for bandwidths
  stopifnot("samplesize has to be greater or equal to the length of the bandwidth collection" = length(samples) >= length(bandwidths))
  stopifnot(is.numeric(bandwidths))
  stopifnot(length(bandwidths) > 0)
  stopifnot(all(bandwidths <= 1) & all(bandwidths >= 1 / length(samples)))
  stopifnot(isTRUE(all(bandwidths > 0)))

  # conditions for subdivisions
  stopifnot(is.numeric(subdivisions))
  stopifnot(length(subdivisions) == 1)
  subdivisions <- ceiling(subdivisions)

  errors <- sapply(bandwidths,
                   cross_validation_error,
                   kernel=kernel,
                   samples=samples,
                   subdivisions=subdivisions)

  bandwidths[which.min(errors)]
}

cross_validation_error <- function(kernel, samples, bandwidth, subdivisions = 1000L) {
  density_estimator <- kernel_density_estimator(kernel, samples, bandwidth, subdivisions)

  squared_l2_norm_estimate <-
    integrate_primitive(
      function(x)
        density_estimator$fun(x) ^ 2,
      lower = density_estimator$support[1],
      upper = density_estimator$support[2],
      subdivisions = subdivisions
    )$value

  num_samples <- length(samples)

  # to use that kernel$fun is vectorised in its argument
  differences <- outer(samples, samples, `-`)
  diag(differences) <- NA_real_
  eval_points <- differences[!is.na(differences)]/bandwidth
  summation <- sum(kernel$fun(eval_points))

  mixed_integral_estimate <- summation/(num_samples*(num_samples-1)*bandwidth)
  squared_l2_norm_estimate - 2*mixed_integral_estimate
}
