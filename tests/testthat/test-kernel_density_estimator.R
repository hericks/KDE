test_that("kernelDensityEstimator returns a valid IntegrableFunction", {
  set.seed(1)
  samples <- rnorm(50)

  kde <- kernelDensityEstimator(rectangular, samples, 0.2, subdivisions=500L)
  expect_error(validate_IntegrableFunction(kde), NA)

  kde <- kernelDensityEstimator(gaussian, samples, 0.05)
  expect_error(validate_IntegrableFunction(kde), NA)
})

test_that("kernelDensityEstimator throws an error for invalid inputs", {
  expect_error(kernelDensityEstimator(sum, 1:10, 1))
  expect_error(kernelDensityEstimator(rectangular, list(1, 2, 3, "a"), 1))
  expect_error(kernelDensityEstimator(rectangular, 1:10, -1))
})

test_that("a KDE evaluates as expected", {
  kde <- kernelDensityEstimator(rectangular, 0, 1)
  eval_points <- c(-1-1e-5, -1, 0, 1, 1+1e5)
  expected_values <- c(0, 0.5, 0.5, 0.5, 0)
  expect_equal(kde$fun(eval_points), expected_values)

  kde <- kernelDensityEstimator(rectangular, c(-0.5, 0.5), 1)
  eval_points <- c(-1.5-1e-5, -1.5, -0.5, 0, 0.5, 1.5, 1.5+1e-5)
  expected_values <- c(0, 0.25, 0.5, 0.5, 0.5, 0.25, 0)
  expect_equal(kde$fun(eval_points), expected_values)
})

test_that("a KDE has the expected support", {
  kde <- kernelDensityEstimator(rectangular, 0, 1)
  expect_equal(kde$support, c(-1,1))

  kde <- kernelDensityEstimator(rectangular, 1, 0.5)
  expect_equal(kde$support, c(0.5,1.5))

  kde <- kernelDensityEstimator(rectangular, -10:10, 0.5)
  expect_equal(kde$support, c(-10.5, 10.5))
})
