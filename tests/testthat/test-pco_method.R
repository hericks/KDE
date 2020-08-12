test_that("kernel must be a valid kernel",{
  not_ker <- function(x){"x"}
  samples <- c(1, 2, 3, 4)
  H_n <- c(0.5, 0.8)
  lambda <- 1
  subdivisions = 10L
  expect_error(pco_method(not_ker,
                          samples = NULL,
                          H_n,
                          lambda,
                          subdivisions))
})

test_that("samples must be numeric with length greater than 0",{
  kernel <- Kernel(epanechnikov_function, c(-1,1), subdivisions=10L)
  H_n <- c(0.5, 0.8)
  lambda <- 1
  subdivisions = 10L
  expect_error(pco_method(kernel,
                          samples = c(1, 2, 3, "a"),
                          H_n,
                          lambda,
                          subdivisions))
  expect_error(pco_method(kernel,
                          samples = NULL,
                          H_n,
                          lambda,
                          subdivisions))
})

test_that("H_n must be numeric and within 1/length(samples) and 1",{
  kernel <- Kernel(epanechnikov_function, c(-1,1), subdivisions=10L)
  samples <- c(1, 2, 3, 4)
  lambda <- 1
  subdivisions = 10L
  #length(H_n) > length(samples)
  expect_error(pco_method(kernel,
                          samples,
                          H_n = c(0.6, 0.7, 0.3, 0.4, 0.5),
                          lambda,
                          subdivisions))
  #contains non-numeric value
  expect_error(pco_method(kernel,
                          samples,
                          H_n = c(0.6, "a"),
                          lambda,
                          subdivisions))
  #length(H_n) = 0
  expect_error(pco_method(kernel,
                          samples,
                          H_n = NULL,
                          lambda,
                          subdivisions))
  #contains value greater than 1
  expect_error(pco_method(kernel,
                          samples,
                          H_n = c(0.6, 1.2),
                          lambda,
                          subdivisions))
  #contains value less than 1/length(samples)
  expect_error(pco_method(kernel,
                          samples,
                          H_n = c(0.1, 0.6),
                          lambda,
                          subdivisions))
  #contains value >= 0
  expect_error(pco_method(kernel,
                          samples,
                          H_n = c(0.6, -0.5),
                          lambda,
                          subdivisions))
})

test_that("lambda has to be numerical scalar",{
  kernel <- Kernel(epanechnikov_function, c(-1,1), subdivisions=10L)
  samples <- c(1, 2, 3, 4)
  H_n <- c(0.5, 0.8)
  subdivisions <- 10L
  expect_error(pco_method(kernel,
                          samples,
                          H_n = NULL,
                          lambda = "a",
                          subdivisions = 100L))
  expect_error(pco_method(kernel,
                          samples,
                          H_n = NULL,
                          lambda = c(1,2),
                          subdivisions = 100L))
  expect_error(pco_method(kernel,
                          samples,
                          H_n = NULL,
                          lambda = TRUE,
                          subdivisions = 100L))
})

test_that("subdivisions must be a single integer",{
  kernel <- Kernel(epanechnikov_function, c(-1,1), subdivisions=10L)
  samples <- c(1, 2, 3, 4)
  H_n <- c(0.5, 0.8)
  lambda <- 1
  expect_error(pco_method(kernel,
                          samples,
                          H_n,
                          lambda = 1,
                          subdivisions = 100))
  expect_error(pco_method(kernel,
                          samples,
                          H_n,
                          lambda = 1,
                          subdivisions = c(1L, 2L)))
  expect_error(pco_method(kernel,
                          samples,
                          H_n,
                          lambda = 1,
                          subdivisions = "a"))
  expect_error(pco_method(kernel,
                          samples,
                          H_n,
                          lambda = 1,
                          subdivisions = FALSE))
})