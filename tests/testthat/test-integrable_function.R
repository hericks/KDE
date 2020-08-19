test_that("the integrable function should be built correctly",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  integ_fun <- IntegrableFunction(rectangular_function,c(-1,1))
  expect_equal(integ_fun$support, c(-1,1))
  expect_true(inherits(integ_fun, "IntegrableFunction"))
  expect_equal(validate_IntegrableFunction(integ_fun), integ_fun)
})

test_that("the ... argument works for integrable functions", {
  expect_equal(IntegrableFunction(dnorm, mean=2)$fun(2), dnorm(0))
  expect_equal(IntegrableFunction(dexp, rate=2)$fun(0), dexp(0, rate=2))
  expect_equal(IntegrableFunction(dunif, max=3)$fun(2.5), 1/3)
})

test_that("the support has to be computed correctly",{
  rectangular_function_a_e_zero <- function(u){
    return(1/2 * 1e5* (abs(u/1e-5) <= 1))
  }
  custom_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  int_fun <- IntegrableFunction(rectangular_function_a_e_zero, NULL)
  int_fun2 <- IntegrableFunction(dnorm, NULL)
  int_fun3 <- IntegrableFunction(dunif, NULL)
  int_fun4 <- IntegrableFunction(dexp, NULL)
  int_fun5 <- IntegrableFunction(custom_den, NULL)
  expect_false(isTRUE(all.equal(1, integrate(rectangular_function_a_e_zero, -1, 1)[[1]])))
  expect_equal(1, integrate(int_fun$fun, int_fun$support[1], int_fun$support[2])[[1]])
  expect_equal(1, integrate(int_fun2$fun, int_fun2$support[1], int_fun2$support[2])[[1]])
  expect_equal(1, integrate(int_fun3$fun, int_fun3$support[1], int_fun3$support[2])[[1]])
  expect_equal(1, integrate(int_fun4$fun, int_fun4$support[1], int_fun4$support[2])[[1]])
  expect_equal(1, integrate(int_fun5$fun, int_fun5$support[1], int_fun5$support[2])[[1]])
})

test_that("the support has to be found or else the object cannot be created",{
  failing_fun <- function(x){rectangular_function <- function(u){
    return(1/2 * (abs(u+1e4+3) <= 1))
  }}
  expect_error(IntegrableFunction(failing_fun, NULL))
})

test_that("the integrable function has to hold a valid support",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  expect_error(IntegrableFunction(rectangular_function, c(0.5, 0.5)))
  expect_error(IntegrableFunction(rectangular_function, c(0.5, -0.5)))
  expect_error(IntegrableFunction(rectangular_function, c("3", "1")))
  expect_error(IntegrableFunction(rectangular_function, 1))
  expect_error(IntegrableFunction(rectangular_function, 1:3))
})

test_that("the integrable function has to be a numeric function",{
  non_numeric_fun <- function(x){"x"}
  expect_error(IntegrableFunction(non_numeric_fun,c(-1,1)))
})

test_that("the function has to be zero outside of support",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  expect_error(IntegrableFunction(rectangular_function,c(-0.5, 0.5)))
})

test_that("non-compact supports or prohibited", {
  expect_error(IntegrableFunction(exp, c(-Inf, Inf)))
})

test_that("validate_IntegrableFunction has to recognize a Integrable function",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  ker <- Kernel(rectangular_function, c(-1,1))
  den <- Density(dnorm, c(-15,15))
  expect_equal(ker, validate_IntegrableFunction(ker))
  expect_equal(den, validate_IntegrableFunction(den))
})

test_that("the subdivisions parameter has to be a positive numeric value",{
  rec_fun <- rectangular$fun
  expect_error(IntegrableFunction(rec_fun, c(-1,1), subdivisions="1.5"))
  expect_error(IntegrableFunction(rec_fun, c(-1,1), subdivisions=c(1,4,5)))
  expect_output(print(IntegrableFunction(rec_fun, c(-1,1), subdivisions=150L)))
})

test_that("the custom printing method works", {
  expected_lines <- "IntegrableFunction\nfunction (x, mean = 0, sd = 1, log = FALSE) .Call(C_dnorm, x, mean, sd, log)\nevaluated at (x, mean = 2), support: [-12,15], subdivisions: 250"
  expect_output(print(IntegrableFunction(dnorm, support=c(-12, 15), subdivisions=250, mean=2)), expected_lines, fixed=TRUE)
})
