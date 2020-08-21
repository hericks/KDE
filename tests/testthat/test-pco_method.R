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
  kernel <- Kernel(epanechnikov_function, c(-1,1))
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
  kernel <- Kernel(epanechnikov_function, c(-1,1))
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
  kernel <- Kernel(epanechnikov_function, c(-1,1))
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

test_that("subdivisions must be a numeric scalar",{
  kernel <- Kernel(epanechnikov_function, c(-1,1))
  samples <- c(1, 2, 3, 4)
  H_n <- c(0.5, 0.8)
  lambda <- 1
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

test_that("pco should return a smaller bandwidth for bigger samplesize", {
  set.seed(50)
  kernel <- epanechnikov

  # Custom density
  f_den_eval <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }

  f_den <- Density(f_den_eval, c(0,1))
  g_den <- Density(dunif, c(0,1))

  # Create sampler from custom density
  custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)

  # Calculate goldenshluger_lepski bandwidth(kernel,
  bandwidth_50 <- cross_validation(kernel, custom_sampler(50), subdivisions = 250L)
  bandwidth_200 <- cross_validation(kernel, custom_sampler(200), subdivisions = 250L)
  bandwidth_500 <- cross_validation(kernel, custom_sampler(500), subdivisions = 250L)
  expect_false(bandwidth_500 > bandwidth_200)
  expect_true(bandwidth_50 > bandwidth_200)
})

