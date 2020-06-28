# Settings
num_samples <- 25
bandwidth <- 0.05
kernel <- gaussian

# Custom density
f_den <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}

# Create sampler from custom density
custom_sampler <- rejection_sampling(f_den, dunif, runif, 2)

# Sample from custom sampler
samples <- custom_sampler(num_samples)

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

