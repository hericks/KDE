test_that("a density should be built correctly",{
  dens_norm <- Density(dnorm, c(-Inf, Inf))
  dens_unif <- Density(dunif, c(-Inf, Inf))
  expect_equal(dens_norm$fun, dnorm)
  expect_equal(dens_unif$fun, dunif)
  expect_equal(dens_norm$support, c(-Inf, Inf))
  expect_equal(dens_unif$support, c(-Inf, Inf))

})
