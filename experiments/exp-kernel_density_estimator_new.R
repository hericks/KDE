
# Custom density
f_den <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}

# Create sampler from custom density
custom_sampler <- rejection_sampling(f_den, dunif, runif, 2)

num_samples <- 100
samples <- custom_sampler(num_samples)
bandwidth <- 0.1

# p_hat <- kernelDensityEstimator(function(x) as.integer(-1 <= x & x <= 1)/2, samples, bandwidth)
p_hat <- kernelDensityEstimator(function(x) as.integer(-1 <= x & x <= 1)*(1-abs(x)), samples, bandwidth)
#p_hat <- kernelDensityEstimator(dnorm, samples, bandwidth)

x <- seq(-3, 3, by=0.005)
y <- p_hat(x)


plot(x, f_den(x), col="red", type="l", ylim=c(0,2.5), xlim=c(-0.5,1.5))
points(samples, integer(length(samples)), pch=".", col="blue")
lines(x, y)
# lines(x, abs(f_den(x) - y))
# integrate(function(x) abs(f_den(x) - p_hat(x)), -0.5, 1.5)
