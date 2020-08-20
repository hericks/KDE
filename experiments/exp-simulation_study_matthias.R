# helper function for plotting
plot_with_confidence_band <- function(x, Y, col) {
  rgb <- col2rgb(col) / 255
  col_alpha <- rgb(rgb[1], rgb[2], rgb[3], 0.2)
  v <- apply(Y, 1, function(x)
    c(mean(x), sd(x)))
  lines(x, v[1,], lwd = 2, col = col)
  polygon(c(x, rev(x)),
          c(v[1,] + v[2,], rev(v[1,] - v[2,])),
          col = col_alpha,
          border = NA)
  lines(x, v[1,] + v[2,], lwd = 1, col = col)
  lines(x, v[1,] - v[2,], lwd = 1, col = col)
}

plot_comparison <- function(density, ...) {
  res <- compare(...)

  x <- res$x[,1]
  y <- res$y

  m <- dim(y)[3]
  del <- max(apply(y, c(1, 3), sd))

  density_values <- density(x)

  # plot results
  par(mar = c(0, 0, 1, 0), ann = FALSE, xaxt = "n", yaxt = "n")

  plot(x, density_values, type = "l", lwd = 2, col = 1, ylim = range(density_values) + del * c(-1, +1))
  grid()

  for (i in 1:m)
    plot_with_confidence_band(x, y[, , i], col = i + 1)
}

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

num_eval_points <- 100
samplers <- list(custom_sampler)
supports <- list(support1=c(-0.5, 5.5))
num_samples <- c(100)
kernels <- list(gaussian, rectangular)
num_bandwidths <- 10L
subdivision_set <- c(1000)
bandwidth_estimators <- list(cross_validation)

library(profvis)

profvis({
  plot_comparison(custom_density_fun,
                  num_eval_points,
                  samplers,
                  supports,
                  num_samples,
                  kernels,
                  num_bandwidths,
                  subdivision_set,
                  bandwidth_estimators,
                  reps=50)
})


