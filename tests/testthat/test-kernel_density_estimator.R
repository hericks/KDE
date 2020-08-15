test_that("KDE returns a valid IntegrableFunction", {
  set.seed(1)
  samples <- rnorm(50)

  kde <- kernel_density_estimator(rectangular, samples, 0.2, subdivisions=500L)
  expect_error(validate_IntegrableFunction(kde), NA)

  kde <- kernel_density_estimator(gaussian, samples, 0.05)
  expect_error(validate_IntegrableFunction(kde), NA)
})

test_that("KDE throws an error for invalid inputs", {
  expect_error(kernel_density_estimator(sum, 1:10, 1))
  expect_error(kernel_density_estimator(rectangular, list(1, 2, 3, "a"), 1))
  expect_error(kernel_density_estimator(rectangular, 1:10, -1))
})

test_that("a KDE evaluates as expected", {
  kde <- kernel_density_estimator(rectangular, 0, 1)
  eval_points <- c(-1-1e-5, -1, 0, 1, 1+1e5)
  expected_values <- c(0, 0.5, 0.5, 0.5, 0)
  expect_equal(kde$fun(eval_points), expected_values)

  kde <- kernel_density_estimator(rectangular, c(-0.5, 0.5), 1)
  eval_points <- c(-1.5-1e-5, -1.5, -0.5, 0, 0.5, 1.5, 1.5+1e-5)
  expected_values <- c(0, 0.25, 0.5, 0.5, 0.5, 0.25, 0)
  expect_equal(kde$fun(eval_points), expected_values)
})

test_that("a KDE has the expected support", {
  kde <- kernel_density_estimator(rectangular, 0, 1)
  expect_equal(kde$support, c(-1,1))

  kde <- kernel_density_estimator(rectangular, 1, 0.5)
  expect_equal(kde$support, c(0.5,1.5))

  kde <- kernel_density_estimator(rectangular, -10:10, 0.5)
  expect_equal(kde$support, c(-10.5, 10.5))
})



test_that("the subdivision parameter has to work properly", {
  # Settings
  num_samples <- 25
  kernel <- rectangular

  # Custom density
  f_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  f_den <- Density(f_den, c(0,1))

  # Create sampler from custom density
  dens_unif <- Density(dunif)
  custom_sampler <- rejection_sampling(f_den, dens_unif, runif, 2)
  samples <- custom_sampler(num_samples)
  expect_output(print(kernel_density_estimator(kernel, samples, 1, 1000L)))
})
