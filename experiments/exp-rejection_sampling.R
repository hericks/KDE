# Custom density
f_den <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}

# Create sampler from custom density
custom_sampler <- rejection_sampling(f_den, dunif, runif, 2)

# Prepare drawing
x <- seq(-0.5, 1.5, by=0.01)
y <- f_den(x)

# Draw density
plot(x, y, type="l", main="Custom density: 1 + sin(2*pi*x)", ylab="density")

# Draw n samples from custom density
n <- 45
points(custom_sampler(n), rep(1, n), col="red", cex=0.5)
