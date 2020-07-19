variance_estim <- function(kernel, samples){
  n <- length(samples)
  l1_kernel <- integrate(function(x){abs(kernel$fun(x))},lower = kernel$support[1], upper = kernel$support[2])
  l2_kernel <- integrate(function(x){kernel$fun(x)^2}, lower = kernel$support[1], upper = kernel$support[2])
  function(h){
    l1_kernel[[1]] * l2_kernel[[1]]/(n*h)
  }
}

#TODO: Faltung und support der Funktion finden (geht das i.A?)
double_kernel_estim <- function(kernel, samples, h, h_ap, grid){
  # center the grid vector
  if(grid %% 2 == 0){
    grid <- grid + 1
  }
  n <- length(samples)
  c <- 1/(h * h_ap * n)
  kernel_h_ap <- kernel_transform(kernel, 0, h_ap)

  kernels_h <- lapply(samples, function(sample) kernel_transform(kernel, sample, h))

  # get the double_kernel_estimator
  fun <- function(vec){
    res_vec <- c()
    for(x in vec){
      res_temp <- c()
      res_temp <- sapply(kernels_h, function(ker_h){
        a <- ker_h$support[1] - abs(kernel_h_ap$support[1])
        b <- ker_h$support[2] + abs(kernel_h_ap$support[2])
        y <- c(seq(from=x, to= a, length.out=(as.integer(grid/2))),
               x,
               seq(from=x, to= b, length.out=as.integer(grid/2)))

        vec <- convolve(kernel_h_ap$fun(y), rev(ker_h$fun(y)), type='open')
        vec <- vec[1:grid]
        vec[as.integer(grid/2) + 1]
      })
      res_vec <- c(res_vec, c * sum(res_temp))
    }
    res_vec
  }

  # get the support of the estimator
  kernels_h_support_lower <- min(sapply(kernels_h, function(x) x$support[1]))
  kernels_h_support_upper <- max(sapply(kernels_h, function(x) x$support[2]))
  #support <- c(kernels_h_support_lower - abs(kernel_h_ap$support[1]),
  #             kernels_h_support_upper + abs(kernel_h_ap$support[2]))


  support <- find_support(fun)
  list("fun"=fun, "support"=support)
}

bias_estim <- function(kernel, samples, h, H_n, kappa = 1.2, var_est, grid_convolve){
  res <- c()
  for (h_ap in H_n) {
    f_double <- double_kernel_estim(kernel, samples, h, h_ap, grid_convolve)
    kde_h_ap <- kernel_density_estimator(kernel, samples, h_ap)

    if (f_double$support[2] < kde_h_ap$support[1] | kde_h_ap$support[2] < f_double$support[1]){
      z1 <- integrate(function(x){f_double(x)^2}, f_double$support[1], f_double$support[2])
      z2 <- integrate(function(x){kde_h_ap(x)^2}, kde_h_ap$support[1], kde_h_ap$support[2])
      l2 <- z1[[1]] + z2[[1]]
    }
    else{
      l2 <- integrate(function(x){(f_double$fun(x) - kde_h_ap$fun(x))^2},
                      lower=min(f_double$support[1], kde_h_ap$support[1]),
                      upper=max(f_double$support[2], kde_h_ap$support[2]))
    }
    res <- c(res, l2[[1]] - kappa * var_est(h_ap))
  }
  max(res)
}

#'
#'
#'@export
goldenshluger_lepski_method <- function(kernel, samples, H_n = NULL, kappa = 1.2, grid_convolve=501L){
  #TODO arg-checks

  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for H_n
  stopifnot("samplesize has to be greater or equal to M"=length(samples) >= length(H_n))
  stopifnot(is.numeric(H_n))
  stopifnot(length(H_n) > 0)
  stopifnot(all(H_n <= 1) & all(H_n >= 1/length(samples)))
  # TODO: Seite 67 Bedingung 3.5 Comte (?), crashes at sample size between 1e7 and 1e8
  stopifnot(!isTRUE(all.equal(1/length(samples), 0)))

  #conditions for kappa
  stopifnot(is.numeric(kappa))
  stopifnot(length(kappa) == 1)
  stopifnot("kappa has to be greater or equal to 1"= kappa >= 1)

  # conditions for grid_convolve
  stopifnot(is.integer(grid_convolve))
  stopifnot(length(grid_convolve) == 1)
  # grid will not work with less than 2 points
  stopifnot("grid_convolve has to be greater or equal to 2"= grid_convolve >= 2)

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
    res_temp <- bias_estim(kernel, samples, h, H_n, kappa, var_est, grid_convolve) + 2 * kappa * var_est(h)
    print(res_temp)
    res <- c(res, res_temp)
  }
  h_est <- H_n[which.min(res)]
}

find_support <- function(fun) {
  testing_points <- c(-10**(10:-10), 10**(-10:10))
  testing_values <- fun(testing_points)

  stopifnot("fun has to return numerical values"=is.numeric(testing_values))
  stopifnot("fun hat to be vectorised"=identical(length(testing_values), length(testing_points)))

  # can't use isFALSE, since all.equal return value is TRUE or a vector of mode "character"
  non_zero_indices <- which(sapply(testing_values, function(x) !isTRUE(all.equal(x, 0))))

  if (length(non_zero_indices) == 0L) return(c(-Inf, Inf))

  lower_index <- min(non_zero_indices) - 1
  upper_index <- max(non_zero_indices) + 1

  lower_bound <- ifelse(lower_index < 1, -Inf, testing_points[lower_index])
  upper_bound <- ifelse(upper_index > length(testing_points), Inf, testing_points[upper_index])

  c(lower_bound, upper_bound)
}
