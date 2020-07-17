variance_estim <- function(kernel, samples){
  n <- length(samples)
  l1_kernel <- integrate(function(x){abs(kernel$fun)},lower = kernel$support[1], upper = kernel$support[2])
  l2_kernel <- integrate(function(x){kernel$fun(x)^2}, lower = kernel$support[1], upper = kernel$support[2])
  function(h){
    l1_kernel*l2_kernel/(n*h)
  }
}

double_kernel_estim <- function(kernel, samples, h, h_ap){
  kernel_h <- kernel_transform(kernel, samples, h)
  kernel_h_ap <- kernel_transform(kernel, 0, h_ap)
  n <- length(samples)
  c <- 1/(h*h_ap * n)
  fun <- function(x){
      #TODO: convolution and support (bias_estim needs support)?
    a <- kernel_h$fun(x)
    b <- kernel_h_ap$fun(x)
    convolve(a,b)

  }
  #TODO
  support <- c(-1,1)
  IntegrableFunction(fun, support)
}

bias_estim <- function(kernel, samples, h, H_n, kappa = 1.2, var_est){
  res <- c()
  for (h_ap in H_n) {
    f_double <- double_kernel_estim(kernel, samples, h, h_ap)
    kde_h_ap <- kernel_density_estimator(kernel, samples, h_ap)
    if (f_double$support[2] < kde_h_ap$support[1] | kde_h_ap$support[2] < f_double$support[1]){
      z1 <- integrate(function(x){f_double(x)^2}, f_double$support[1], f_double$support[2])
      z2 <- integrate(function(x){kde_h_ap(x)^2}, kde_h_ap$support[1], kde_h_ap$support[2])
      l2 <- z1+z2
    }
    else{
      l2 <- integrate(function(x){(f_double(x) - kde_h_ap(x))^2},
                      lower=min(f_double$support[1], kde_h_ap$support[1]),
                      upper=max(f_double$support[2], kde_h_ap$support[2]))
    }
    res <- c(res, l2 -kappa*var_est(h_ap) )
  }
  max(res)
}

goldenshluger_lepski_method <- function(kernel, samples, H_n = NULL, kappa = 1.2){
  #TODO arg-checks

  # Kernel conditions
  stopifnot("the kernel has to be valid" = validate_Kernel(kernel))

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for H_n
  stopifnot("samplesize has to be greater or equal to M"=length(samples) >= length(H_n))
  stopifnot(is.numeric(H_n))
  stopifnot(length(H_n) > 0)
  stopifnot(all(H_n <= 1) & all(H_n >= 1/length(samples)))
  #TODO Seite 67 3.5 Compte

  #conditions for kappa
  stopifnot("kappa has to be greater or equal to 1"= kappa >= 1)
  stopifnot(is.numeric(kappa))
  stopifnot(length(kappa)== 1)

  if(is.null(H_n)) {
    H_n <- c()
    n <-length(samples)
    for (m in (1:as.integer(n))){
      H_n <- c(H_n, m/n)
    }
  }

  res <- c()
  var_est <- variance_estim(kernel, samples)
  for (h in H_n){
    res_temp <- bias_estim(kernel, samples, h, M, kappa, var_est(h)) + 2*kappa*var_est(h)
    res <- c(res, res_temp)
  }
  h_est <- which.min(res)/M
}



