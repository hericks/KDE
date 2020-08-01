test_that("a KDE is valid IntegrableFunction", {
  kde <- kernelDensityEstimator(rectangular, 1:10, 0.2)

  expect_error(validate_IntegrableFunction(kde), NA)
})
