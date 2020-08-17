#' Penalized Comparison to Overfitting
#'
#' @description The PCO method is used to estimate an optimal bandwidth for
#'   kernel density estimation from a given set of bandwidths.
#'
#' @param kernel S3 object of class \code{\link{Kernel}}; the kernel to use for
#'   the estimator
#' @param samples numeric vector; the observations.
#' @param bandwidths strictly positive numeric vector; the bandwidth set from
#'   which the bandwidth with the least estimated risk will be selected.
#' @param lambda positive numeric scalar; a tuning parameter.
#' @param subdivisions positive numeric scalar; subdivisions parameter
#'   internally passed to \code{\link{integrate_primitive}}.
#'
#' @return The estimated optimal bandwidth contained in the bandwidth set.
#'
#' @details The PCO method aims to minimize an upper bound for the mean
#'   integrated squared error (MISE) of a kernel density estimator. The MISE is
#'   defined as the expectation of the squared L2-Norm of the difference between
#'   estimator and (unknown) true density.
#'
#'   \code{pco_method} internally uses a criterion function to calculate the PCO
#'   criterion value, approximating the risk. Subsequently the bandwidth with
#'   the minimal criterion value is selected.
#'
#'   The popular bias-/variance-decomposition is used. The bias term still
#'   depends on the unknown density. Thus, a comparison of the estimator with an
#'   associated bandwidth to the overfitting one, namely the estimator with the
#'   smallest bandwidth, is used to estimate the bias term itself.
#'
#'   Further a penalty term is computed as the sum of two variances,
#'   particularly the variance of the risk decomposition and the variance of the
#'   bias term estimation. During the calculation the tuning parameter
#'   \code{lambda} is used. The recommended value for \code{lambda} is 1.
#'
#'   The PCO criterion is given by the sum of the comparison to overfitting and
#'   the penalty term, thus the procedure tries to find a balance between those
#'   terms. Therefore, it is comprehensible why this method is called penalized
#'   comparison to overfitting.
#'
#'   For more information see the linked papers below.
#'
#' @seealso \code{\link{kernel_density_estimator}} for more information about
#'   kernel density estimators, \code{\link{cross_validation}} and
#'   \code{\link{goldenshluger_lepski}} for more automatic bandwidth-selection
#'   algorithms.
#'
#' @source \href{https://arxiv.org/abs/1607.05091v2}{Estimator selection: a new
#'   method with applications to KDEs}, Lacour \[2017\]
#' @source \href{https://arxiv.org/abs/1902.01075}{Numerical performance of
#'   PCO for multivariate KDEs}, Varet \[2019\]
#'
#' @include kernel.R
#' @include kernel_density_estimator.R
#' @include logarithmic_bandwidth_set.R
#'
#' @export
pco_method <- function(kernel,
                       samples,
                       bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10),
                       lambda = 1,
                       subdivisions = 100L) {

  # Argchecks are being made in pco_crit, hence the pco_crit function can be used without pco_method, but not the other way around.
  res_tuple <- pco_crit(kernel, samples, bandwidths, lambda, subdivisions)
  h_est <- bandwidths[which.min(res_tuple$risk)]
}

pco_crit <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), lambda = 1, subdivisions = 100L) {
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

  # conditions for lambda
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)

  # conditions for subdivisions
  stopifnot(is.numeric(subdivisions))
  stopifnot(length(subdivisions) == 1)
  subdivisions <- ceiling(subdivisions)

  res <- c()
  h_min <- min(bandwidths)
  f_h_min_est <- kernel_density_estimator(kernel, samples, h_min, subdivisions)

  for (h in bandwidths) {
    f_h_est <- kernel_density_estimator(kernel, samples, h, subdivisions=subdivisions)

    bias_estim <-
      integrate_primitive(
        function(x) {
          (f_h_min_est$fun(x) - f_h_est$fun(x)) ^ 2
        },
        lower = min(f_h_min_est$support[1], f_h_est$support[1]),
        upper = max(f_h_min_est$support[2], f_h_est$support[2]), subdivisions=subdivisions
      )$value

    pen_function <- penalty_term(kernel, samples, h_min, h, lambda)
    l_pco <- bias_estim + pen_function
    res <- c(res, l_pco)
  }
  list(bandwidth_set = bandwidths, risk = res)
}

penalty_term <- function(kernel, samples, h_min, h, lambda, subdivisions = 1000L) {
  n <- length(samples)
  ker_h_min <- kernel_transform(kernel, 0, h_min, subdivisions)
  ker_h <- kernel_transform(kernel, 0, h, subdivisions)

  # calculate the variance term
  l2_h_kernel <-
    integrate_primitive(function(x) {
      ker_h$fun(x) ^ 2
    },
    lower = ker_h$support[1],
    upper = ker_h$support[2], subdivisions=subdivisions)$value
  h_var_term <- lambda * l2_h_kernel / n


  z <- integrate_primitive(
    function(x) {
      (ker_h_min$fun(x) - ker_h$fun(x)) ^ 2
    },
    lower = min(ker_h_min$support[1], ker_h$support[1]),
    upper = max(ker_h_min$support[2], ker_h$support[2]),
    subdivisions = subdivisions
  )$value

  bias_term_pre <- z
  bias_term <- bias_term_pre / n

  penalty <- h_var_term - bias_term
}
