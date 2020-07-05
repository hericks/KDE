test_that("kernel needs to be a kernel function",{
  non_kernel_fun <- function(x) 1/x
  expect_error(kernelTransform(non_kernel_fun,0,1))
})

test_that("bandwidth has to be positive numerical scalar",{
  non_numerical <- "a"
  expect_error(kernelTransform(gaussian,0,non_numerical))
  neg_bandwidth <- -1
  expect_error(kernelTransform(gaussian,0,neg_bandwidth))
  non_scalar <- 1:2
  expect_error(kernelTransform(gaussian,0,non_scalar))
})

test_that("observation has to be numerical scalar",{
  non_numerical_obs <- "a"
  expect_error(kernel_transform(gaussian,non_numerical_obs,1))
  non_scalar_obs <- 1:2
  expect_error(kernelTransform(gaussian,non_scalar_obs,1))
})


test_that("bandwidth is working",{
  h1 <- 2
  h2 <- 0.5
  obs <- 0
  triangular_stretched <- kernelTransform(triangular, obs, h1)
  triangular_compressed <- kernelTransform(triangular, obs, h2)
  expect_true(triangular_stretched(0) < triangular(0))
  expect_true(triangular_compressed(0) > triangular(0))
})

test_that("observation is shifting the kernel",{
  h <- 1
  obs <- 2
  triangular_shifted <- kernelTransform(triangular, obs, h)
  expect_equal(triangular(0), triangular_shifted(2))
  expect_true(triangular(2) < triangular_shifted(2))
  expect_true(triangular(0) > triangular_shifted(0))
})
