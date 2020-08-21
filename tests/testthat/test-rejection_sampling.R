# "example for which algorithm was designed"
test_that("the samples are distributed according to the given cdf with density function 1+sin(2*pi*x)",{
  custom_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  f_den <- Density(custom_den, c(0,1))
  g_den <- Density(dunif)

  custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)
  y <- custom_sampler(1e6)
  expect_gte(min(y), 0)
  expect_lte(max(y), 1.0)
  expect_gte(max(y), 0.9)
  expect_lte(min(y), 0.1)
})

# test of a "trivial" example where you can use g_den = f_den and M = 1
test_that("if den is density for uniform distribution, sampels are uniformily distributed",{
  den <- Density(dunif)
  custom_sampler <- rejection_sampling(den, den, runif, 1)
  y <- custom_sampler(1e6)
  expect_gte(min(y), 0)
  expect_lte(max(y), 1.0)
  expect_gte(max(y), 0.9)
  expect_lte(min(y), 0.1)
  expect_lt(abs(mean(y) - 0.5), 0.1)
  expect_lt(var(y), 1/6)
  expect_gt(var(y), 1/24)
})

test_that("if den is density for normal distribution, sampels are normally distributed",{
  den <- Density(dnorm)
  custom_sampler <- rejection_sampling(den, den, rnorm, 1)
  y <- custom_sampler(1e6)
  expect_gte(min(y), den$support[1])
  expect_lte(max(y), den$support[2])
  expect_lt(abs(mean(y)), 0.1)
  eps <- 0.2
  expect_lt(var(y), 1+eps)
  expect_gt(var(y), 1-eps)
})

test_that("f_den has to be a Density object",{
  den <- Density(dunif)
  expect_error(rejection_sampling(dunif, den, runif, 1))
  expect_error(rejection_sampling(den, dunif, runif, 1))
  expect_error(rejection_sampling("dunif", den, runif, 1))
  expect_error(rejection_sampling(den, "dunif",  runif, 1))
})

test_that("f_den hat to be less then or equal M * g_den",{
  f_den <- Density(dunif)
  g_den <- Density(dnorm)
  expect_error(rejection_sampling(f_den, g_den, rnorm, 1))
  expect_error(rejection_sampling(g_den, f_den, runif, 1))
})

test_that("g has to be a numerical function",{
  den <- Density(dunif)
  expect_error(rejection_sampling(den, den, "runif", 1))
  expect_error(rejection_sampling(den, den, 1:3, 1))
})


test_that("M has to be a numerical scalar",{
  den <- Density(dunif)
  expect_error(rejection_sampling(den, den, runif, "1"))
  expect_error(rejection_sampling(den, den, runif, 1:3))
})

