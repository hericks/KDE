#' Goldenshluger-Lepski Method
#'
#' @description The Goldenshluger-Lepski Method is used to estimate a optimal bandwidth for kernel density estimation from a given set of bandwidths.
#'
#' @param kernel kernel function as an S3 object of the class \link[KDE:Kernel]{Kernel}.
#' @param samples A numerical vector of observations.
#' @param bandwidths The bandwidth set from which the bandwidth with the least risk will be selected. The \code{goldenshluger_lepski} function will try to set up a suitable bandwidth set using \code{\link[KDE:logarthmic_bandwidth_set]{logarithmic_bandwidth_set}} if \code{NULL} is passed.
#' @param kappa A tuning parameter. It has to be a numerical value with length 1. The minimal admissible value is 1. The recommendation is to set kappa = 1.2.
#' @param subdivisions A integer vector of length 1 used for the subdivisions parameter of the builtin R-function \code{\link[stats:integrate]{integrate}}. The default value is set to 100L.
#'
#' @return A numerical vector of length 1 containing the bandwidth with minimal risk.
#'
#' @details This method is an implementation of the Goldenshluger-Lepski method for bandwith estimation. The aim is to minimize the risk of a kernel density estimator (KDE).
#' The risk function is given by the expected value of the integrated square error of the KDE with a given bandwith and the desired density that is matching the distibution of the samples (which we are trying to estimate with the KDE). \cr
#' By applying this method the risk will be decomposed into a bias and a variance term, where the
#' bias term needs to be estimated, because it includes a dependency of the real function that we are
#' trying to estimate with the KDE. The variance term is simply a bound for the variance of the KDE (created with a
#' bandwith \code{h}), which is tuned by a parameter \code{kappa}.
#' The bias term will be estimated by using a double kernel approach.
#' The method then selects the bandwidth with the minimal associated risk. \cr
#' For further information about the Goldenshluger-Lepski method for bandwidth estimation, see "Nonparametric Estimation" by Fabienne Comte.
#'
#' @seealso
#' \itemize{\code{\link[KDE:pco_method]{PCO Method}},
#' \code{\link[KDE:cross_validation]{Cross-Validation Method}} to see other methods for estimating bandwidths.}
#' \itemize{\code{\link[KDE:Kernel]{Kernel}} to see the definiotion of a kernel.}
#' \itemize{\code{\link[KDE:kernel_density_estimator]{Kernel Density Estimator}} to see the functionality of the Kernel density estimation.}
#'
#' @source
#' \itemize{\href{https://spartacus-idh.com/030.html}{Nonparametric Estimation}, Comte `[`2017`]`, ISBN: 978-2-36693-030-6}
#'
#' @include kernel.R
#' @include kernel_density_estimator.R
#' @include logarithmic_bandwidth_set.R
#'
#' @export
goldenshluger_lepski <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), kappa = 1.2, subdivisions = 100L) {
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
  stopifnot(isTRUE(all(bandwidths > 0)))

  # conditions for kappa
  stopifnot(is.numeric(kappa))
  stopifnot(length(kappa) == 1)
  stopifnot("kappa has to be greater or equal to 1" = kappa >= 1)

  # conditions for subdivisions
  stopifnot(is.integer(subdivisions))
  stopifnot(length(subdivisions) == 1)

  # Used for calculation of V(h)
  squared_l1_norm_kernel <- integrate(function(x) abs(kernel$fun(x)), lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)[[1]]**2
  squared_l2_norm_kernel <- integrate(function(x) kernel$fun(x)^2, lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)[[1]]

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

    squared_l2_norm <- integrate(function(x) {
        (f_h_h2(x) - kde_h2$fun(x))^2
      },
      lower = min(f_h_h2_support, kde_h2$support),
      upper = max(f_h_h2_support, kde_h2$support),
      subdivisions = subdivisions)[[1]]

    V_h2 <- squared_l1_norm_kernel*squared_l2_norm_kernel/(length(samples)*h2)

    ret <- max(ret, squared_l2_norm - kappa*V_h2)
  }

  ret
}

# Calculate convolution of f with g, where both function have a support inside [x_min, x_max]
conv <- function(f, g, x_min, x_max, N = 300) {
  x <- seq(x_min, x_max, len = N)
  s <- (x_max - x_min) / N
  y_out <- s * convolve(f(x), g(rev(x)), conj=TRUE, type="open")
  x_out <- seq(2*x_min+s, 2*x_max-s, len = 2*N-1) # offset s is artifact of discretization
  list(x = x_out, y = y_out)
}
