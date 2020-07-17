# Settings
num_samples <- 25
bandwidth <- 0.05
kernel <- gaussian

# Custom density
f_den_eval <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}

f_den <- Density(f_den_eval, c(0,1))
g_den <- Density(dunif, c(0,1))

# Create sampler from custom density
custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)

# Sample from custom sampler
samples <- custom_sampler(num_samples)

# Create KDE
p_hat <- kernelDensityEstimator(kernel, samples, bandwidth)$fun

# Plot the density in red
x <- seq(-3, 3, by=0.005)
y <- p_hat(x)
plot(x, f_den_eval(x), col="red", type="l", ylim=c(0,2.5), xlim=c(-0.5,1.5))

# Plot the samples in blue
points(samples, integer(length(samples)), pch=".", col="blue")

# Plot the estimator in black
lines(x, y)

