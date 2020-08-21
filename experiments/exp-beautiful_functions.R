custom_density_fun <- function(x) {
  res <- double(length(x))
  res[0 <= x & x <= 1] <- x[0 <= x & x <= 1]
  res[1 <= x & x <= 2] <- 1
  res[2 <= x & x <= 3] <- x[2 <= x & x <= 3] - 1
  res[3 <= x & x <= 5] <- 5 - x[3 <= x & x <= 5]
  res/5
}

custom_density <- Density(custom_density_fun, c(0, 5))
unif_density <- Density(function(x) { dunif(x, 0, 5) }, c(0, 5))

custom_sampler <- rejection_sampling(custom_density,
                                     unif_density,
                                     function(n) {
                                       runif(n, 0, 5)
                                     }, 2)
