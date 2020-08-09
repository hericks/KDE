cross_validation_error <- function(kernel, samples, bandwidth, subdivisions=100L) {
  density_estimator <- kernelDensityEstimator(kernel, samples, bandwidth, subdivisions)

  # TODO: Usable error if integration fails (increase number of subdivisions)
  squared_l2_norm_estimate <-
    integrate(
      function(x)
        density_estimator$fun(x) ^ 2,
      lower = density_estimator$support[1],
      upper = density_estimator$support[2],
      subdivisions = subdivisions
    )[[1]]

  num_samples <- length(samples)

  # to use that kernel$fun is vectorised in its argument
  differences <- outer(samples, samples, `-`)
  diag(differences) <- NA_real_
  eval_points <- differences[!is.na(differences)]/bandwidth
  summation <- sum(kernel$fun(eval_points))

  mixed_integral_estimate <- summation/(num_samples*(num_samples-1)*bandwidth)
  squared_l2_norm_estimate - 2*mixed_integral_estimate
}
