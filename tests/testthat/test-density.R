test_that("a density should be built correctly",{
  dens_norm <- Density(dnorm, c(-15, 15))
  dens_unif <- Density(dunif, c(0, 1))
  custom_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  dens_custom <- Density(custom_den, c(0,1))

  expect_equal(dens_norm$support, c(-15, 15))
  expect_equal(dens_unif$support, c(0, 1))
  expect_equal(dens_custom$support, c(0,1))
  expect_true(inherits(dens_norm, "Density"))
  expect_true(inherits(dens_unif, "Density"))
  expect_true(inherits(dens_custom, "Density"))
})

test_that("the ... argument works for densities", {
  expect_equal(Density(dnorm, mean=2)$fun(2), dnorm(0))
  expect_equal(Density(dunif, min=-1, max=1)$fun(0), 1/2)
})


test_that("density has to hold a valid support", {
  expect_error(Density(dnorm, c(0.5,0.5)))
  expect_error(Density(dnorm, c(0.5,-0.5)))
  expect_error(Density(dnorm, c(-0.5,0.5)))
  expect_error(Density(dnorm, c("3", "1")))
  expect_error(Density(dnorm, 1))
  expect_error(Density(dnorm, 1:3))
})

test_that("the integral of a density over the real numbers equals one",{
  div_integral_fun <- function(x) 1
  finite_integral_fun <- function(x) 1*(x <= 100)
  non_integrable_fun <- function(x) 1/x

  expect_error(Density(div_integral_fun,c(0,2)))
  expect_error(Density(finite_integral_fun,c(0,100)))
  expect_error(Density(non_integrable_fun,c(0,1)))
})

test_that("density function has to be a numeric function",{
  non_numeric_fun <- function(x){"x"}

  expect_error(Density(non_numeric_fun,c(-1,1)))
})


test_that("density has to be non-negative",{
  neg_valued_function <- function(x){-1}
  neg_valued_function2 <- function(x){
    sapply(x, function(x){
      if((x >= -1) & (x < -0.5)) return(-1)
      else{
          return(1 * (abs(x) <= 1))
      }})
  }
  neg_valued_function3 <- function(x){
    sapply(x, function(x){
      if((x >= -0.9) & (x < -0.4)) return(-1)
      else{
        return(1 * (abs(x) <= 1))
      }})
  }
  expect_error(Density(neg_valued_function, c(-Inf, Inf)))
  expect_error(Density(neg_valued_function2, c(-1,1)))
  expect_error(Density(neg_valued_function3, c(-1,1)))
})

test_that("the custom printing method works", {
  expected_lines <- "Density\nfunction (x, mean = 0, sd = 1, log = FALSE) .Call(C_dnorm, x, mean, sd, log)\nevaluated at (x, mean = 2), support: [-12,15], subdivisions: 250"
  expect_output(print(Density(dnorm, support=c(-12, 15), subdivisions=250, mean=2)), expected_lines, fixed=TRUE)
})
