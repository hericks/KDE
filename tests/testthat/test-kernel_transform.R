test_that("kernel needs to be a kernel function",{
  non_kernel_fun <- function(x) 1/x
  non_kernel_fun <- class(c("Kernel", "IntegrableFunction"))
  expect_error(kernel_transform(non_kernel_fun,0,1))
})

test_that("kernel_transform does not change the kernel object",{
  gaussian_temp <- gaussian
  kernel_transform(gaussian_temp, 1, 5)
  expect_true(identical(gaussian_temp, gaussian))
})

test_that("kernel_transform returns valid kernel object", {
  expect_error(validate_Kernel(kernel_transform(gaussian, 1, 2)), NA)
  expect_error(validate_Kernel(kernel_transform(rectangular, 10, 0.2)), NA)
})

test_that("bandwidth has to be positive numerical scalar",{
  non_numerical <- "a"
  expect_error(kernel_transform(gaussian,0,non_numerical))
  neg_bandwidth <- -1
  expect_error(kernel_transform(gaussian,0,neg_bandwidth))
  non_scalar <- 1:2
  expect_error(kernel_transform(gaussian,0,non_scalar))
})

test_that("sample has to be numerical scalar",{
  non_numerical_sample <- "a"
  expect_error(kernel_transform(gaussian,non_numerical_sample,1))
  non_scalar_sample <- 1:2
  expect_error(kernel_transform(gaussian,non_scalar_sample,1))
})


test_that("bandwidth is working",{
  h1 <- 2
  h2 <- 0.5
  sample <- 0
  triangular_stretched <- kernel_transform(triangular, sample, h1)
  triangular_compressed <- kernel_transform(triangular, sample, h2)
  expect_true(triangular_stretched$fun(0) < triangular$fun(0))
  expect_true(triangular_compressed$fun(0) > triangular$fun(0))
})

test_that("sample is shifting the kernel",{
  h <- 1
  sample <- 2
  triangular_shifted <- kernel_transform(triangular, sample, h)
  expect_equal(triangular$fun(0), triangular_shifted$fun(2))
  expect_true(triangular$fun(2) < triangular_shifted$fun(2))
  expect_true(triangular$fun(0) > triangular_shifted$fun(0))
})
