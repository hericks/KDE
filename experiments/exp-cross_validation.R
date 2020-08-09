# Settings
num_samples <- 100
bandwidths <- c(1, 0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01)
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

# Calculate cross-validation errors
errors <- sapply(bandwidths, function(h) cross_validation_error(kernel, samples, h, 100))

# Print optimal bandwidth
paste("Selected bandwidth:", bandwidths[which.min(errors)])

# Create KDEs
estimators <- lapply(bandwidths, function(h) kernelDensityEstimator(kernel, samples, h, 100)$fun)

# Plot the density in red
x <- seq(-3, 3, by=0.005)
plot(x, f_den_eval(x), col="red", type="l", ylim=c(0,2.5), xlim=c(-0.5,1.5))

# Plot the samples in blue
points(samples, integer(length(samples)), pch=".", col="blue")

# Plot best estimator
lines(x, estimators[[which.min(errors)]](x))

