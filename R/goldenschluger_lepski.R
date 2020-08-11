goldenschluger_lepski <- function(kernel, samples, bandwidths, subdivisions = 100L) {
  # Used for calculation of V(h)
  squared_l1_norm_kernel <- integrate(function(x) abs(kernel$fun(x)), lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)[[1]]**2
  squared_l2_norm_kernel <- integrate(function(x) kernel$fun(x)^2, lower=kernel$support[1], upper=kernel$support[2], subdivisions = subdivisions)[[1]]

  penalties <- numeric(0)
  for (h in bandwidths) {
    A <- A(kernel, samples, bandwidths, h, subdivisions, squared_l1_norm_kernel, squared_l2_norm_kernel)
    V <- squared_l1_norm_kernel*squared_l2_norm_kernel/(length(samples)*h)
    penalties[length(penalties) + 1] <- A + 2*V
  }

  bandwidths[which.min(penalties)]
}

A <- function(kernel, samples, bandwidths, h, subdivisions = 100L, squared_l1_norm_kernel, squared_l2_norm_kernel) {
  # Stretched kernel
  k_h <- kernel_transform(kernel, 0, h)

  # Calculate content of supremum for each bandwidth h2
  ret <- -Inf
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

    ret <- max(ret, squared_l2_norm - V_h2)
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
