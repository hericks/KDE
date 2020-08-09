variance_estim <- function(kernel, samples){
  n <- length(samples)
  l1_kernel <- integrate(function(x){abs(kernel$fun(x))},lower = kernel$support[1], upper = kernel$support[2])
  l2_kernel <- integrate(function(x){kernel$fun(x)^2}, lower = kernel$support[1], upper = kernel$support[2])
  function(h){
    l1_kernel[[1]]^2 * l2_kernel[[1]]/(n*h)
  }
}

conv <- function(f, g, x_min, x_max, N = 300) {
  x <- seq(x_min, x_max, len = N)
  s <- (x_max - x_min) / N
  y_out <- s * convolve(f(x), g(rev(x)), conj=TRUE, type="open")
  x_out <- seq(2*x_min+s, 2*x_max-s, len = 2*N-1) # offset s is artifact of discretization
  list(x = x_out, y = y_out)
}

double_kernel_estim <- function(kernel, samples, h, h_ap, grid_size){
  n <- length(samples)
  kernel_h_ap <- kernel_transform(kernel, 0, h_ap)
  kernel_h <- kernel_density_estimator(kernel, samples, h)

  x_min <- kernel_h$support[1] - abs(kernel_h_ap$support[1])
  x_max <- kernel_h$support[2] + abs(kernel_h_ap$support[2])

  convolution <- conv(kernel_h_ap$fun, kernel_h$fun, x_min, x_max, grid_size)
  fun <- approxfun(convolution$x, convolution$y)

  support <- c(x_min, x_max)
  list("fun"=fun, "support"=support)
}

bias_estim <- function(kernel, samples, h, H_n, kappa = 1.2, var_est, grid_size){
  res <- c()
  for (h_ap in H_n) {
    f_double <- double_kernel_estim(kernel, samples, h, h_ap, grid_size)
    kde_h_ap <- kernel_density_estimator(kernel, samples, h_ap)

    if (f_double$support[2] < kde_h_ap$support[1] | kde_h_ap$support[2] < f_double$support[1]){
      z1 <- integrate(function(x){f_double(x)^2}, f_double$support[1], f_double$support[2])
      z2 <- integrate(function(x){kde_h_ap(x)^2}, kde_h_ap$support[1], kde_h_ap$support[2])
      l2 <- z1[[1]] + z2[[1]]
    }
    else{
      l2 <- integrate(function(x){(f_double$fun(x) - kde_h_ap$fun(x))^2},
                      lower=min(f_double$support[1], kde_h_ap$support[1]),
                      upper=max(f_double$support[2], kde_h_ap$support[2]),
                      subdivisions=2000)
      l2 <- l2[[1]]
    }

    print(l2 - kappa * var_est(h_ap))
    res <- c(res, l2 - kappa * var_est(h_ap))
  }
  res[res < 0] <- 0
  max(res)
}

#'
#'
#'@export
goldenshluger_lepski_method_2 <- function(kernel, samples, H_n = NULL, kappa = 1.2, grid_size=501L){

  if(is.null(H_n)) {
    H_n <- c()
    n <-length(samples)
    for (m in (1:as.integer(n))){
      H_n <- c(H_n, m/n)
    }
  }

  #TODO arg-checks

  # conditions for kernel
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # conditions for samples
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # conditions for H_n
  stopifnot("samplesize has to be greater or equal to M"= length(samples) >= length(H_n))
  stopifnot(is.numeric(H_n))
  stopifnot(length(H_n) > 0)
  stopifnot(all(H_n <= 1) & all(H_n >= 1/length(samples)))
  #stopifnot(!isTRUE(all.equal(1/length(samples), 0)))
  stopifnot(isTRUE(all(H_n > 0)))

  # conditions for kappa
  stopifnot(is.numeric(kappa))
  stopifnot(length(kappa) == 1)
  stopifnot("kappa has to be greater or equal to 1"= kappa >= 1)

  # conditions for grid_size
  stopifnot(is.integer(grid_size))
  stopifnot(length(grid_size) == 1)
  # grid will not work with less than 2 points
  stopifnot("grid_convolve has to be greater or equal to 2"= grid_size >= 2)



  res <- c()
  var_est <- variance_estim(kernel, samples)
  for(h in H_n){
    print("-------------------------")
    res_temp <- bias_estim(kernel, samples, h, H_n, kappa, var_est, grid_size) + 2 * kappa * var_est(h)
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
