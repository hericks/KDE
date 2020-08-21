raw_kde <- function(kernel, samples, bandwidth=1) {
  kernel_eval <- kernel$fun

  estimator_eval <- function(x) {
    stopifnot("input has to be numeric"=is.numeric(x))

    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel_eval((x - x0)/bandwidth)/bandwidth
    }
    ret/length(samples)
  }

  support <- bandwidth*kernel$support + range(samples)

  list(
    fun=estimator_eval,
    support=support
  )
}

# Settings
fixed = 1000
subdivisions = 1000L

# Estimators (multiple sample sets)
bandwidth <- 1

dnorm_sample_sets <- list(
  rnorm(10), rnorm(50), rnorm(100), rnorm(250), rnorm(500), rnorm(1000)
)

estimators <- lapply(
  dnorm_sample_sets, function(samples) {
    raw_kde(epanechnikov, samples, bandwidth)
  }
)

# Estimators (multiple bandwidths)
bandwidths <- c(1, 0.5, 0.1, 0.05, 0.025, 0.01, 0.005)

samples <- rnorm(1000)

estimators <- lapply(
  bandwidths, function(h) {
    raw_kde(triangular, samples, h)
  }
)

# Error Analysis
primitive_abs_errors <- double(length(estimators))
for (i in seq_along(estimators)) {
  estimator <- estimators[[i]]
  I <- integrate_primitive(integrand = estimator$fun,
                           lower = estimator$support[1],
                           upper = estimator$support[2],
                           subdivisions = fixed)$value

  primitive_abs_errors[i] <- abs(I - 1)
}

base_abs_errors <- double(length(estimators))
for (i in seq_along(estimators)) {
  print("OK")
  estimator <- estimators[[i]]
  I <- integrate(f = estimator$fun,
                 lower = estimator$support[1],
                 upper = estimator$support[2],
                 subdivisions = subdivisions)$value

  base_abs_errors[i] <- abs(I - 1)
}

summary(primitive_abs_errors)
summary(base_abs_errors)

# Performance Analysis
primitive_integration <- function() {
  for (estimator in estimators) {
    integrate_primitive(integrand = estimator$fun,
                        lower = estimator$support[1],
                        upper = estimator$support[2],
                        subdivisions = fixed)$value
  }
}

base_integration <- function() {
  for (estimator in estimators) {
    integrate(f = estimator$fun,
              lower = estimator$support[1],
              upper = estimator$support[2],
              subdivisions = subdivisions)
  }
}

microbenchmark::microbenchmark(primitive_integration(), base_integration(), times=1L)
