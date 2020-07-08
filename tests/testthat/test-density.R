test_that("a density should be built correctly",{
  dens_norm <- Density(dnorm, c(-Inf, Inf))
  dens_unif <- Density(dunif, c(-Inf, Inf))
  custom_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  dens_custom <- Density(custom_den, c(0,1))

  expect_equal(dens_norm$fun, dnorm)
  expect_equal(dens_unif$fun, dunif)
  expect_equal(dens_custom$fun, custom_den)
  expect_equal(dens_norm$support, c(-Inf, Inf))
  expect_equal(dens_unif$support, c(-Inf, Inf))
  expect_equal(dens_custom$support, c(0,1))
  expect_true(inherits(dens_norm, "Density"))
  expect_true(inherits(dens_unif, "Density"))
  expect_true(inherits(dens_custom, "Density"))
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

  expect_error(Density(div_integral_fun,c(-Inf,Inf)))
  expect_error(Density(finite_integral_fun,c(-Inf,Inf)))
  expect_error(Density(non_integrable_fun,c(-Inf,Inf)))
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
}
)
