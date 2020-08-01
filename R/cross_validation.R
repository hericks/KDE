cross_validation_error <- function(kernel, samples, bandwidth, subdivisions=NULL) {
  if (is.null(subdivisions)) {
    subdivisions = 35*sqrt(length(samples))/sqrt(bandwidth)
    print(subdivisions)
  }

  density_estimator <- kernelDensityEstimator(kernel, samples, bandwidth)

  squared_l2_norm_estimate <-
    integrate(
      function(x)
        density_estimator$fun(x) ^ 2,
      lower = density_estimator$support[1],
      upper = density_estimator$support[2],
      subdivisions = subdivisions
    )[[1]]

  num_samples <- length(samples)

  temp <- 0
  for(i in seq_along(samples)) {
    for(j in seq_along(samples)) {
      if (i == j) next

      temp <- temp + kernel$fun((samples[i] - samples[j])/bandwidth)
    }
  }

  mixed_integral_estimate <- temp/(num_samples*(num_samples-1)*bandwidth)

  squared_l2_norm_estimate - 2*mixed_integral_estimate
}
