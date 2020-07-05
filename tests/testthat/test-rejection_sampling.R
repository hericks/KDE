test_that("the samples are distributed according to the given cdf with density function 1+sin(2*pi*x)",{
  f_den <- function(x) {
    ret <- 1 + sin(2*pi*x)
    ret[x < 0 | 1 < x] <- 0
    ret
  }
  custom_sampler <- rejection_sampling(f_den, dunif, runif, 2)
  y <- custom_sampler(1e6)
  expect_gte(min(y), 0)
  expect_lte(max(y), 1.0)
  expect_gte(max(y), 0.9)
  expect_lte(min(y), 0.1)
})
