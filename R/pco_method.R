penalty_term <- function(kernel, samples, h_min, h, lambda, subdivisions = 100L) {
  n <- length(samples)
  ker_h_min <- kernel_transform(kernel, 0, h_min, subdivisions)
  ker_h <- kernel_transform(kernel, 0, h, subdivisions)

  # calculate the variance term
  l2_h_kernel <-
    integrate(function(x) {
      ker_h$fun(x) ^ 2
    },
    lower = ker_h$support[1],
    upper = ker_h$support[2])
  h_var_term <- lambda * l2_h_kernel[[1]] / n

  if (ker_h_min$support[2] < ker_h$support[1] |
      ker_h$support[2] < ker_h_min$support[1]) {
    z1 <-
      integrate(function(x) {
        ker_h_min$fun(x) ^ 2
      }, ker_h_min$support[1], ker_h_min$support[2])
    z2 <-
      integrate(function(x) {
        ker_h$fun(x) ^ 2
      }, ker_h$support[1], ker_h$support[2])
    bias_term_pre <- z1[[1]] + z2[[1]]
  }
  else{
    z <- integrate(
      function(x) {
        (ker_h_min$fun(x) - ker_h$fun(x)) ^ 2
      },
      lower = min(ker_h_min$support[1], ker_h$support[1]),
      upper = max(ker_h_min$support[2], ker_h$support[2])
    )
    bias_term_pre <- z[[1]]
  }
  bias_term <- bias_term_pre / n

  penalty <- h_var_term - bias_term
}


#' Penalized Comparison To Overfitting criterion calculation
#' @description pco_crit calculates the associated risk for each bandwidth in a given bandwidth set.
#'
#' @param kernel kernel function as an S3 object of the class \link[KDE:Kernel]{Kernel}.
#' @param samples A numerical vector of observations to base the construction of the estimator.
#' @param H_n The bandwidth set from which the bandwidth with the least risk according to the PCO criterion will be derived. The pco_methdod function will try to set up a suitable bandwidth set if \code{NULL} is passed.
#' @param lambda A tuning parameter. It has to be a numerical value with length 1. The criterion blows up for lambda < 0, therefore the optimal value for lambda is a positiv real number. The recommendation is to set lambda = 1.
#' @param subdivisions A integer vector of length 1 used for the subdivisions parameter of the builtin R-function \code{\link[stats:integrate]{integrate}}.
#'
#' @details For more information about the Penalized Comparison to Overfitting method, see [pco_method].
#'
#' @return A list containing the bandwidth set and associated risk values according to the PCO criterion.
#'
#' @seealso \code{\link[KDE:pco_method]{pco_method}} for more information about the PCO method.
#'
#'#' @include kernel.R
#'
#' @export
pco_crit <- function(kernel, samples, H_n, lambda, subdivisions = 100L) {
  # conditions for kernel
  tryCatch({
    validate_Kernel(kernel)
  }, error = "the kernel has to be valid")

  # conditions for samples
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for H_n
  stopifnot("samplesize has to be greater or equal to M" = length(samples) >= length(H_n))
  stopifnot(is.numeric(H_n))
  stopifnot(length(H_n) > 0)
  stopifnot(all(H_n <= 1) & all(H_n >= 1 / length(samples)))
  #stopifnot(!isTRUE(all.equal(1/length(samples), 0)))
  stopifnot(isTRUE(all(H_n > 0)))


  # conditions for lambda
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)

  # conditions for subdivisions
  stopifnot(is.integer(subdivisions))
  stopifnot(length(subdivisions) == 1)

  res <- c()
  h_min <- min(H_n)
  f_h_min_est <- kernel_density_estimator(kernel, samples, h_min, subdivisions)

  for (h in H_n) {
    f_h_est <- kernel_density_estimator(kernel, samples, h, subdivisions)
    if (f_h_min_est$support[2] < f_h_est$support[1] |
        f_h_est$support[2] < f_h_min_est$support[1]) {
      z1 <-
        integrate(function(x) {
          f_h_min_est$fun(x) ^ 2
        }, f_h_min_est$support[1], f_h_min_est$support[2], subdivisions=subdivisions)
      z2 <-
        integrate(function(x) {
          f_h_est$fun(x) ^ 2
        }, f_h_est$support[1], f_h_est$support[2], subdivisions=subdivisions)
      bias_estim <- z1[[1]] + z2[[1]]
    }
    else{
      bias_estim <-
        integrate(
          function(x) {
            (f_h_min_est$fun(x) - f_h_est$fun(x)) ^ 2
          },
          lower = min(f_h_min_est$support[1], f_h_est$support[1]),
          upper = max(f_h_min_est$support[2], f_h_est$support[2]), subdivisions=subdivisions
        )
      bias_estim <- bias_estim[[1]]
    }
    pen_function <- penalty_term(kernel, samples, h_min, h, lambda)
    l_pco <- bias_estim + pen_function
    res <- c(res, l_pco)
  }
  list(bandwidth_set = H_n, risk = res)
}


#' Penalized Comparison to Overfitting
#'
#' @description The PCO method is used to estimate a optimal bandwidth for kernel density estimation from a given set of bandwidths.
#'
#' @param kernel kernel function as an S3 object of the class \link[KDE:Kernel]{Kernel}.
#' @param samples A numerical vector of observations.
#' @param H_n The bandwidth set from which the bandwidth with the least risk according to the PCO criterion will be derived. The \code{pco_method} function will try to set up a suitable bandwidth set if \code{NULL} is passed.
#' @param lambda A tuning parameter. It has to be a numerical value with length 1. The criterion blows up for lambda < 0, therefore the optimal value for lambda is a positiv real number. The recommendation is to set lambda = 1.
#' @param subdivisions A integer vector of length 1 used for the subdivisions parameter of the builtin R-function \code{\link[stats:integrate]{integrate}}.
#'
#' @return A numerical vector of length 1 containing the bandwidth with minimal risk according to the PCO criterion, given by \code{\link[pco_crit]{pco_crit}}.
#'
#' @details pco_method uses \code{\link[KDE:pco_crit]{pco_crit}} to calculate the PCO criterion value, approximating the risk.
#' For each bandwidth given in \code{H_n}, the method selects the bandwidth with the minimal associated risk.
#'
#' The risk function is given by the expected value of the integrated square error of the KDE with
#' a given bandwith and the desired density that is matching the distibution of the samples (which we
#'                                                                                           are trying to estimate with the KDE).
#' By applying this method the risk will be decomposed into a bias and a variance term, where the
#' bias term needs to be estimated, because it includes a dependency of the real function that we are
#' trying to estimate with the KDE. The variance term is simply a bound for the variance of the KDE (created with a
#' bandwith \code{h}), which is tuned by a parameter \code{lambda}.\cr
#' The bias term will be estimated by a comparison of the KDE of a bandwidth \code{h} to the KDE of the
#' smallest bandwidth \code{h_min} out of the given bandwidth set \code{H_n}. \cr
#' A penalty term consisting of the sum of the two variance terms is introduced, including the variance from the risk decomposition and the one from
#' the bias term estimation.
#' The PCO criterion is given by the sum of the comparison to overfitting and the penalty term.
#' Now it is comprehensible why this method is called Penalized Comparison to Overfitting. \cr
#' In the end the method compares the KDE of each bandwidth in \code{H_n} to the KDE of the smallest bandwidth out
#' of the given bandwidth collection. Finally the bandwidth with minimal associated risk will be returned. \cr
#' For more information see the linked papers below.
#'
#'
#'
#'
#' @seealso
#' \itemize{\code{\link[KDE:pco_crit]{pco_crit}} for getting the associated criterion values.}
#' \itemize{\code{\link[KDE:goldenshluger_lepski_method]{Goldenshluger-Lepski method}},
#' \code{\link[KDE:cross_validation]{Cross-Validation method}} to see other methods for estimating bandwidths.}
#' \itemize{\code{\link[KDE:Kernel]{Kernel}} to see the definiotion of a kernel.}
#' \itemize{\code{\link[KDE:kernel_density_estimator]{Kernel Density Estimator}} to see the functionality of the Kernel density estimation.}
#'
#' @source \href{https://arxiv.org/abs/1607.05091v2}{Lacour `[`2017`]`}
#' @source \href{https://arxiv.org/abs/1902.01075}{Varet `[`2019`]`}
#'
#' @include kernel.R
#' @export
pco_method <- function(kernel,
                       samples,
                       H_n = NULL,
                       lambda = 1,
                       subdivisions = 100L) {
  if (is.null(H_n)) {
    num_samples <- length(samples)
    H_n <- log(1 - seq(1, 1/num_samples, length.out=20))/log(1 - 1/num_samples)
    H_n <- H_n[is.finite(H_n)] - min(H_n)
    H_n <- 1/num_samples + (1 - 1/num_samples)*H_n/max(H_n)
  }

  # Argchecks are being made in pco_crit, hence the pco_crit function can be used without pco_method, but not the other way around.

  res_tuple <- pco_crit(kernel, samples, H_n, lambda, subdivisions)
  h_est <- H_n[which.min(res_tuple$risk)]
}
