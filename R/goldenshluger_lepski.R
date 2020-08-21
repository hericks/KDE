#' Goldenshluger-Lepski Method
#'
#' @description The Goldenshluger-Lepski method is used to estimate an optimal
#'   bandwidth for kernel density estimation from a given set of bandwidths.
#'
#' @param kernel S3 object of class \code{\link{Kernel}}; the kernel to use for
#'   the estimator
#' @param samples numeric vector; the observations.
#' @param bandwidths strictly positive numeric vector; the bandwidth set from
#'   which the bandwidth with the least estimated risk will be selected.
#' @param kappa numeric scalar greater 1; a tuning parameter.
#' @param subdivisions positive numeric scalar; subdivisions parameter
#'   internally passed to \code{\link{integrate_primitive}}.
#'
#' @details The Goldenshluger-Lepski method aims to minimize an upper bound for
#'   the mean integrated squared error (MISE) of a kernel density estimator. The
#'   MISE is defined as the expectation of the squared L2-Norm of the difference
#'   between estimator and (unknown) true density.
#'
#'   This methods works with the popular bias-/variance-decomposition. A
#'   double-kernel approach is used for an estimator of the bias term as it
#'   still depends on the unknown density.
#'
#'   The estimator used for an upper bound of the variance depends on the tuning
#'   parameter \code{kappa}. The recommended value for \code{kappa} is 1.2.
#'
#'   Subsequently the bandwidth with the minimal associated risk is selected.
#'
#' @return The estimated optimal bandwidth contained in the bandwidth set.
#'
#' @seealso \code{\link{kernel_density_estimator}} for more information about
#'   kernel density estimators, \code{\link{pco_method}} and
#'   \code{\link{cross_validation}} for more automatic bandwidth-selection
#'   algorithms.
#'
#' @source \href{https://spartacus-idh.com/030.html}{Nonparametric Estimation},
#'   Comte \[2017\], ISBN: 978-2-36693-030-6
#'
#' @include kernel.R
#' @include kernel_density_estimator.R
#' @include logarithmic_bandwidth_set.R
#'
#' @importFrom stats approxfun
#'
#' @export
goldenshluger_lepski <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), kappa = 1.2, subdivisions = 1000L) {
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

  # conditions for kappa
  stopifnot(is.numeric(kappa))
  stopifnot(length(kappa) == 1)
  stopifnot("kappa has to be greater or equal to 1" = kappa >= 1)

  # conditions for subdivisions
  stopifnot(is.numeric(subdivisions))
  stopifnot(length(subdivisions) == 1)
  subdivisions <- ceiling(subdivisions)

  # Used for calculation of V(h)
  squared_l1_norm_kernel <- integrate_primitive(function(x) abs(kernel$fun(x)), lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)$value**2
  squared_l2_norm_kernel <- integrate_primitive(function(x) kernel$fun(x)^2, lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)$value

  penalties <- numeric(0)
  for (h in bandwidths) {
    A <- A(kernel, samples, bandwidths, h, kappa, subdivisions, squared_l1_norm_kernel, squared_l2_norm_kernel)
    V <- squared_l1_norm_kernel*squared_l2_norm_kernel/(length(samples)*h)
    penalties[length(penalties) + 1] <- A + 2*kappa*V
  }

  bandwidths[which.min(penalties)]
}

A <- function(kernel, samples, bandwidths, h, kappa, subdivisions = 100L, squared_l1_norm_kernel, squared_l2_norm_kernel) {
  # Stretched kernel
  k_h <- kernel_transform(kernel, 0, h)

  # Calculate content of supremum for each bandwidth h2
  ret <- 0
  for(h2 in bandwidths) {
    kde_h2 <- kernel_density_estimator(kernel, samples, h2, subdivisions = subdivisions)
    k_h2 <- kernel_transform(kernel, 0, h2)

    joint_support <- range(k_h$support, k_h2$support)
    discrete_conv <- conv(k_h$fun, k_h2$fun, joint_support[1], joint_support[2])
    continuous_conv <- approxfun(discrete_conv$x, discrete_conv$y)

    f_h_h2 <- function(x) {
      res <- 0
      for (sample in samples) {
        temp_res <- continuous_conv(x - sample)
        temp_res[is.na(temp_res)] <- 0
        res <- res + temp_res
      }
      res/length(samples)
    }

    f_h_h2_support <- range(discrete_conv$x) + range(samples)

    squared_l2_norm <- integrate_primitive(function(x) {
        (f_h_h2(x) - kde_h2$fun(x))^2
      },
      lower = min(f_h_h2_support, kde_h2$support),
      upper = max(f_h_h2_support, kde_h2$support),
      subdivisions = subdivisions)$value

    V_h2 <- squared_l1_norm_kernel*squared_l2_norm_kernel/(length(samples)*h2)

    ret <- max(ret, squared_l2_norm - kappa*V_h2)
  }

  ret
}

# Calculate convolution of f with g, where both function have a support inside [x_min, x_max]
#' @importFrom stats convolve
conv <- function(f, g, x_min, x_max, N = 300) {
  x <- seq(x_min, x_max, len = N)
  s <- (x_max - x_min) / N
  y_out <- s * convolve(f(x), g(rev(x)), conj=TRUE, type="open")
  x_out <- seq(2*x_min+s, 2*x_max-s, len = 2*N-1) # offset s is artifact of discretization
  list(x = x_out, y = y_out)
}
