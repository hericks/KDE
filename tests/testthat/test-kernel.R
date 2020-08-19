test_that("a kernel should be built correctly",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  ker <- Kernel(rectangular_function,c(-1,1))
  expect_equal(ker$support, c(-1,1))
  expect_true(inherits(ker, "Kernel"))

  # find_borders has to work on rectangular_function
  ker <- Kernel(rectangular_function, NULL)
  # based on how find_borders works, it will find c(-10,10) as borders
  expect_equal(ker$support, c(-10,10))
  expect_true(inherits(ker, "Kernel"))

  # the find_borders has to work on small bandwidth
  rectangular_function_a_e_zero <- function(u){
    return(1/2 * 1000 * (abs(u/0.001) <= 1))
  }
  ker <- Kernel(rectangular_function_a_e_zero,NULL)
  expect_equal(0, integrate(rectangular_function_a_e_zero, -1, 1)[[1]])
  expect_equal(ker, Kernel(rectangular_function_a_e_zero,c(-0.01,0.01)))
  expect_true(inherits(ker, "Kernel"))
})

test_that("the ... argument works for kernels", {
  expect_equal(Kernel(dnorm, mean=2)$fun(2), dnorm(0))
  expect_equal(Kernel(dunif, min=-1, max=1)$fun(0), 1/2)
})

test_that("kernel has to hold a valid support", {
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  expect_error(Kernel(rectangular_function, c(0.5,0.5)))
  expect_error(Kernel(rectangular_function, c(0.5,-0.5)))
  expect_error(Kernel(rectangular_function, c(-0.5,0.5)))
  expect_error(Kernel(rectangular_function, c("3", "1")))
  expect_error(Kernel(rectangular_function, 1))
  expect_error(Kernel(rectangular_function, 1:3))
})

test_that("the integral of a kernel over the real numbers equals one",{
  div_integral_fun <- function(x) 1
  finite_integral_fun <- function(x) 1*(x <= 100)
  non_integrable_fun <- function(x) 1/x

  expect_error(Kernel(div_integral_fun,c(0,2)))
  expect_error(Kernel(finite_integral_fun,c(0,100)))
  expect_error(Kernel(non_integrable_fun,c(0,1)))
})

test_that("kernel function has to be a numeric function",{
  non_numeric_fun <- function(x){"x"}

  expect_error(Kernel(non_numeric_fun,c(-1,1)))
})

test_that("the custom printing method works", {
  expected_output <- "Kernel\nfunction (x, mean = 0, sd = 1, log = FALSE) .Call(C_dnorm, x, mean, sd, log)\nevaluated at (x, mean = 2), support: [-12,15], subdivisions: 250"
  expect_output(print(Kernel(dnorm, support=c(-12, 15), subdivisions=250, mean=2)), expected_output, fixed=TRUE)
})
