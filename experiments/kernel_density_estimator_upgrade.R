kernel_density_estimator_update <- function(kernel, samples, bandwidth=1, subdivisions=1000L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  # Subdivisions conditions
  subdivisions <- subdivisions
  stopifnot("Entry 'subdivisions' must be numeric"=is.numeric(subdivisions))
  stopifnot("Entry 'subdivisions' must be positive"=subdivisions > 0)

  kernel_eval <- kernel$fun

  formal_estimator_eval <- function(x) {
    stopifnot("input has to be numeric"=is.numeric(x))

    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel_eval((x - x0)/bandwidth)/bandwidth
    }
    ret/length(samples)
  }

  support <- bandwidth*kernel$support + range(samples)

  found <- FALSE
  num_points <- 500
  while(!found) {
    x <- seq(support[1], support[2], length.out = num_points)
    approx_eval <- approxfun(x, formal_estimator_eval(x))

    print("Appriximated")

    practical_eval <- function(x) {
      res <- approx_eval(x)
      res[is.na(res)] <- 0
      res
    }

    current_integral <- integrate_primitive(practical_eval, support[1], support[2], subdivisions = 10*num_points)$value
    print(current_integral)
    if (abs(current_integral - 1) < 0.01)
      break

    num_points <- 2*num_points
  }

  IntegrableFunction(practical_eval, support, subdivisions=subdivisions)
}

