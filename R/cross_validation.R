#' @export
cross_validation <- function(kernel, samples, bandwidths = logarithmic_bandwidth_set(1/length(samples), 1, 10), subdivisions = 100L) {
  if (is.null(bandwidths)) {
    num_samples <- length(samples)
    bandwidths <- log(1 - seq(1, 1/num_samples, length.out=20))/log(1 - 1/num_samples)
    bandwidths <- bandwidths[is.finite(bandwidths)] - min(bandwidths)
    bandwidths <- 1/num_samples + (1 - 1/num_samples)*bandwidths/max(bandwidths)
  }
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
  #stopifnot(!isTRUE(all.equal(1/length(samples), 0)))
  stopifnot(isTRUE(all(bandwidths > 0)))

  # conditions for subdivisions
  stopifnot(is.integer(subdivisions))
  stopifnot(length(subdivisions) == 1)

  errors <- sapply(bandwidths, function(h) cross_validation_error(kernel, samples, h, subdivisions = subdivisions))
  bandwidths[which.min(errors)]
}

cross_validation_error <- function(kernel, samples, bandwidth, subdivisions = 100L) {

  density_estimator <- kernel_density_estimator(kernel, samples, bandwidth, subdivisions)

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
