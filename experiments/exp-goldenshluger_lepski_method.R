library(KDE)
# Settings
num_samples <- 25
bandwidth_set <- 1/(1:num_samples)
kernel <- rectangular

# Custom density
f_den <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}
f_den <- Density(f_den, c(0,1))

# Create sampler from custom density
dens_unif <- Density(dunif)
custom_sampler <- rejection_sampling(f_den, dens_unif, runif, 2)

# Sample from custom sampler
samples <- custom_sampler(num_samples)

# goldenshluger_lepski bandwidth estimation
bandwidth <- goldenshluger_lepski_method(kernel, samples, bandwidth_set, grid_convolve=2L)
print(bandwidth)
# Create KDE
p_hat <- kernelDensityEstimator(kernel, samples, bandwidth)

# Plot the density in red
x <- seq(-3, 3, by=0.005)
y <- p_hat(x)
plot(x, f_den(x), col="red", type="l", ylim=c(0,2.5), xlim=c(-0.5,1.5))

# Plot the samples in blue
points(samples, integer(length(samples)), pch=".", col="blue")

# Plot the estimator in black
lines(x, y)
