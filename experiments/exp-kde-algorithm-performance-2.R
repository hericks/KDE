library(profvis)

samples <- rnorm(100)
bandwidth <- 1/length(samples)

profvis({
  estimator1 <- kernel_density_estimator(gaussian, samples, bandwidth, subdivisions = 1000L)
  estimator2 <- kernel_density_estimator_update(gaussian, samples, bandwidth, subdivisions = 1000L)

  integrate_primitive(estimator1$fun, estimator1$support[1], estimator1$support[2], subdivisions = 1000L)$value
  integrate_primitive(estimator2$fun, estimator2$support[1], estimator2$support[2], subdivisions = 1000L)$value
})

profvis({
  goldenshluger_lepski(gaussian, samples)
})

samples <- rnorm(10000)
bandwidth <- 0.0001

profvis({
  estimator <- kernel_density_estimator(rectangular, samples, bandwidth, subdivisions = 1000L)
  alternative_estimator <- kernel_density_estimator_alternative(rectangular, samples, bandwidth, subdivisions = 1000L)

  integrate_primitive(estimator$fun, estimator$support[1], estimator$support[2], subdivisions = 5000L)
  integrate_primitive(alternative_estimator$fun, alternative_estimator$support[1], alternative_estimator$support[2], subdivisions = 5000L)
})

microbenchmark::microbenchmark(integrate_primitive(estimator$fun, estimator$support[1], estimator$support[2]),
                               integrate_primitive(alternative_estimator$fun, alternative_estimator$support[1], alternative_estimator$support[2]))

goldenshluger_lepski(gaussian, rnorm(10000))
