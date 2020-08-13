# Settings
num_samples <- 1000
bandwidth_set <- c(0.5, 0.001)
kernel <- gaussian
bandwidth_set_2 <- c()
for (m in (1:as.integer(40))){
  bandwidth_set_2 <- c(bandwidth_set_2, m/40)
}

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

# pco bandwidth estimation
bandwidth <- goldenshluger_lepski(kernel, samples, bandwidth_set, subdivisions=250L)
print(bandwidth)

# Create KDE
p_hat <- kernel_density_estimator(kernel, samples, bandwidth, subdivisions=250L)

# Plot the density in red
x <- seq(-3, 3, by=0.005)
y <- p_hat$fun(x)
plot(x, f_den$fun(x), col="red", type="l", ylim=c(0,2.5), xlim=c(-0.5,1.5))

# Plot the samples in blue
points(samples, integer(length(samples)), pch=".", col="blue")

# Plot the estimator in black
lines(x, y)

