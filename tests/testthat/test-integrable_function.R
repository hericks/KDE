test_that("the integrable function should be built correctly",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  integ_fun <- IntegrableFunction(rectangular_function,c(-1,1))
  expect_equal(integ_fun$fun, rectangular_function)
  expect_equal(integ_fun$support, c(-1,1))
  expect_true(inherits(integ_fun, "IntegrableFunction"))
  expect_equal(validate_IntegrableFunction(integ_fun), integ_fun)
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

test_that("The integral of the absolute function must integrate to a finite value",{
  diverging_fun <- function(x) 1
  diverging_fun2 <- function(x) x
  expect_error(IntegrableFunction(diverging_fun, c(-Inf, Inf)))
  expect_error(IntegrableFunction(diverging_fun2, c(-Inf, Inf)))
})

test_that("the integrable function has to be integrable over its support",{
  non_integrable <- function(x) 1/x
  expect_error(IntegrableFunction(non_integrable, c(-Inf,Inf)))
})

test_that("validate_IntegrableFunction has to recognize a Integrable function",{
  rectangular_function <- function(u){
    return(1/2 * (abs(u) <= 1))
  }
  ker <- Kernel(rectangular_function, c(-1,1))
  den <- Density(dnorm, c(-Inf,Inf))
  expect_equal(ker, validate_IntegrableFunction(ker))
  expect_equal(den, validate_IntegrableFunction(den))
})