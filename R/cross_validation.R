#' Cross Validation
#'
#' @description The Cross Validation method is used to estimate a optimal bandwidth for kernel density estimation from a given set of bandwidths.
#'
#' @param kernel kernel function as an S3 object of the class \link[KDE:Kernel]{Kernel}.
#' @param samples A numerical vector of observations.
#' @param bandwidths The bandwidth set from which the bandwidth with the least estimated risk will be selected. The \code{cross_validation} function will try to set up a suitable bandwidth set using \code{\link[KDE:logarthmic_bandwidth_set]{logarithmic_bandwidth_set}} if \code{NULL} is passed.
#' @param subdivisions A integer vector of length 1 used for the subdivisions parameter of the builtin R-function \code{\link[stats:integrate]{integrate}}.
#'
#' @return A numerical vector of length 1 containing the bandwidth with minimal cross validation error.
#'
#' @details This method is an implementation of the cross validation method for bandwith estimation. The aim is to minimize the risk of a kernel density estimator (KDE). \cr
#' For each bandwidth \code{h} given in \code{bandwidths}, \code{cross_validation} approximates the term of the risk that is dependent of the KDE.
#' The risk function is given by the expected value of the integrated square error of the KDE with a given bandwith and the desired density that is matching the distibution of the samples (which we are trying to estimate with the KDE).
#' The method then selects the bandwidth with the minimal associated risk. \cr
#' For more information about cross validation for bandwidth estimation, see "Nonparametric Estimation" by Fabienne Comte.
#'
#' @seealso
#' \itemize{\code{\link[KDE:goldenshluger_lepski]{Goldenshluger-Lepski method}},
#' \code{\link[KDE:pco_method]{ PCO method}} to see other methods for estimating bandwidths.}
#' \itemize{\code{\link[KDE:Kernel]{Kernel}} to see the definiotion of a kernel.}
#' \itemize{\code{\link[KDE:kernel_density_estimator]{Kernel Density Estimator}} to see the functionality of the Kernel density estimation.}
#'
#'@source
#' \itemize{\href{https://spartacus-idh.com/030.html}{Nonparametric Estimation}, Comte `[`2017`]`, ISBN: 978-2-36693-030-6}
#'
#' @include kernel.R
#' @include kernel_density_estimator.R
#' @include logarithmic_bandwidth_set.R
#' @export
cross_validation <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), subdivisions = 1000L) {
  # conditions for kernel
  tryCatch({
    validate_Kernel(kernel)
  }, error = "the kernel has to be valid")

  # conditions for samples
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for H_n
  stopifnot("samplesize has to be greater or equal to the length of the bandwidth collection" = length(samples) >= length(bandwidths))
  stopifnot(is.numeric(bandwidths))
  stopifnot(length(bandwidths) > 0)
  stopifnot(all(bandwidths <= 1) & all(bandwidths >= 1 / length(samples)))
  #stopifnot(!isTRUE(all.equal(1/length(samples), 0)))
  stopifnot(isTRUE(all(bandwidths > 0)))

  # conditions for subdivisions
  stopifnot(is.integer(subdivisions))
  stopifnot(length(subdivisions) == 1)

  errors <- sapply(bandwidths, function(h) cross_validation_error(kernel, samples, h, subdivisions = subdivisions))
  bandwidths[which.min(errors)]
}

cross_validation_error <- function(kernel, samples, bandwidth, subdivisions = 1000L) {

  density_estimator <- kernel_density_estimator(kernel, samples, bandwidth, subdivisions)

  # TODO: Usable error if integration fails (increase number of subdivisions)
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
