kernel_density_estimator_alternative <- function(kernel, samples, bandwidth=1, subdivisions=1000L) {
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

  x <- seq(support[1], support[2], length.out = subdivisions)
  formal_values <- formal_estimator_eval(x)

  approximation <- approxfun(x, formal_values)

  practical_eval <- function(x) {
    res <- approximation(x)
    res[is.na(res)] <- 0
    res
  }

  IntegrableFunction(practical_eval, support, subdivisions=20*subdivisions)
}
