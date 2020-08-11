#' Penalized-Comparison-to-Overfitting-Method (PCO)
#'
#' @description The penalty term is essential to estimate the variance, it will be used in
#'
#' @param kernel The (vectorised) kernel function to use for the construction of
#'   the estimator satisfying [is_kernel].
#' @param samples A numerical vector to base the construction of the estimator
#' @param h_min The minimum of the bandwidths grid. It must be a numerical vector with length one.
#' @param h The bandwith we want to be evaluated for the penalty term (?). It also must be a numerical vetor, with lenghth one.
#' @param lambda A tuning parameter. It has to be a numerical value with length 1. The risk blows up for lambda < 0, therefor the optimal value for lambda is a positiv real number. The recommendation is to set lambda=1.
#'
#' @export
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

pco_run <- function(kernel, samples, H_n, lambda, subdivisions = 100L) {
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
        }, f_h_min_est$support[1], f_h_min_est$support[2])
      z2 <-
        integrate(function(x) {
          f_h_est$fun(x) ^ 2
        }, f_h_est$support[1], f_h_est$support[2])
      bias_estim <- z1[[1]] + z2[[1]]
    }
    else{
      bias_estim <-
        integrate(
          function(x) {
            (f_h_min_est$fun(x) - f_h_est$fun(x)) ^ 2
          },
          lower = min(f_h_min_est$support[1], f_h_est$support[1]),
          upper = max(f_h_min_est$support[2], f_h_est$support[2])
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
#'@export
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


  res_tuple <- pco_run(kernel, samples, H_n, lambda, subdivisions)
  h_est <- H_n[which.min(res_tuple$risk)]
}
