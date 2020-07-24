#' @export
penalty_term <- function(kernel, samples, h_min, h, lambda){
  n <- length(samples)
  ker_h_min <- kernel_transform(kernel, 0, h_min)
  ker_h <- kernel_transform(kernel, 0, h)

  # calculate the variance term
  l2_h_kernel <- integrate(function(x){ker_h$fun(x)^2}, lower = ker_h$support[1], upper = ker_h$support[2])
  h_var_term <- lambda * l2_h_kernel[[1]] / n

  if (ker_h_min$support[2] < ker_h$support[1] | ker_h$support[2] < ker_h_min$support[1]){
    z1 <- integrate(function(x){ker_h_min$fun(x)^2}, ker_h_min$support[1], ker_h_min$support[2])
    z2 <- integrate(function(x){ker_h$fun(x)^2}, ker_h$support[1], ker_h$support[2])
    bias_term_pre <- z1[[1]] + z2[[1]]
  }
  else{
    z <- integrate(function(x){(ker_h_min$fun(x) - ker_h$fun(x))^2},
                            lower=min(ker_h_min$support[1], ker_h$support[1]),
                            upper=max(ker_h_min$support[2], ker_h$support[2]))
    bias_term_pre <- z[[1]]
  }
  bias_term <- bias_term_pre / n

  penalty <- h_var_term - bias_term
}

#' Penalized Comparison to Overfitting
#'
#'@export
pco_method <- function(kernel, samples, H_n = NULL, lambda = 1){
  if(is.null(H_n)) {
    H_n <- c()
    n <-length(samples)
    for (m in (1:as.integer(n))){
      H_n <- c(H_n, m/n)
    }
  }

  #TODO arg-checks

  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # TODO: conditions for H_n
  stopifnot("samplesize has to be greater or equal to M"=length(samples) >= length(H_n))
  stopifnot(is.numeric(H_n))
  stopifnot(length(H_n) > 0)
  stopifnot(all(H_n <= 1) & all(H_n >= 1/length(samples)))
  # TODO: Seite 67 Bedingung 3.5 Comte (?), crashes at sample size between 1e7 and 1e8
  stopifnot(!isTRUE(all.equal(1/length(samples), 0)))

  # conditions for lambda
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)



  res <- c()
  h_min <- min(H_n)
  f_h_min_est <- kernel_density_estimator(kernel, samples, h_min)

  for(h in H_n){
    f_h_est <- kernel_density_estimator(kernel, samples, h)
    if (f_h_min_est$support[2] < f_h_est$support[1] | f_h_est$support[2] < f_h_min_est$support[1]){
      z1 <- integrate(function(x){f_h_min_est$fun(x)^2}, f_h_min_est$support[1], f_h_min_est$support[2])
      z2 <- integrate(function(x){f_h_est$fun(x)^2}, f_h_est$support[1], f_h_est$support[2])
      bias_estim <- z1[[1]] + z2[[1]]
    }
    else{
      bias_estim <- integrate(function(x){(f_h_min_est$fun(x) - f_h_est$fun(x))^2},
                      lower=min(f_h_min_est$support[1], f_h_est$support[1]),
                      upper=max(f_h_min_est$support[2], f_h_est$support[2]))
      bias_estim <- bias_estim[[1]]
    }
    pen_function <- penalty_term(kernel, samples, h_min, h, lambda)
    l_pco <- bias_estim + pen_function
    res <- c(res, l_pco)
  }
  h_est <- H_n[which.min(res)]
}
