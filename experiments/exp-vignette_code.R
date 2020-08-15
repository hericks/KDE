# besser funktionierenden Kern nehmen
# 3 plots mit rectangular  verschiedene bandweiten/sample sizes (1x sehr wenig samples (4), 1x 25 samples, 1x 150 samples)
# bandbreitenwahl nochmal 3 plots mit verschiedenen bandbreiten (1x over 1x under und 1x relativ gut)
# rest schreiben

# Custom density
f_dens <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}
support_density <- c(0,1)
# constructing the Density object
dens <- Density(f_dens, support_density, subdivisions=100L)
print(dens)


# or build a Kernel object yourself: this is the sigmoid function
f_ker <- function(x){
  return(2/pi * (1/(exp(x) + exp(-x))))
}
# A function is sufficiently small outside of a compact space, such that the numeric integration outside of this space will be zero.
# We need to chose a finite Interval as support for some bandwidth estimators to work properly.
support_sigmoid <- c(-20,20)
# constructing the Kernel object
sigmoid_kernel <- Kernel(f_ker, support_sigmoid, subdivisions=100L)
print(sigmoid_kernel)


# using a builtin Kernel
gaussian_kernel <- gaussian
print(gaussian_kernel)

# Create sampler from custom density
g_den <- Density(dunif, c(0,1))
custom_sampler <- rejection_sampling(dens, g_den, runif, 2)
# Sample from custom sampler
num_samples <- 100
samples <- custom_sampler(num_samples)

# setting up bandwidth sets
bandwidth_set <- logarithmic_bandwidth_set(from=1/length(samples), to=1, length.out=10)

# bandwidth estimation
# cross-validation method
cv_gaussian <- cross_validation(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
cv_sigmoid <- cross_validation(sigmoid_kernel, samples, bandwidths=bandwidth_set, subdivisions=300L)
cat("bandwidth for gaussian kernel: ", cv_gaussian, "\nbandwidth for sigmoid kernel: ", cv_sigmoid)

# goldenshluger-lepski method
gl_gaussian <- goldenshluger_lepski(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
gl_sigmoid <- goldenshluger_lepski(sigmoid_kernel, samples, bandwidths=bandwidth_set, subdivisions=300L)
cat("bandwidth for gaussian kernel: ", gl_gaussian, "\nbandwidth for sigmoid kernel: ", gl_sigmoid)

# Penalized Comparison to Overfitting
pco_gaussian <- pco_method(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
pco_sigmoid <- pco_method(sigmoid_kernel, samples, bandwidths=bandwidth_set, subdivisions=300L)
cat("bandwidth for gaussian kernel: ", pco_gaussian, "\nbandwidth for sigmoid kernel: ", pco_sigmoid)


# Kernel Density Estimation


kde_cv_gaussian <- kernel_density_estimator(gaussian_kernel, samples, cv_gaussian, subdivisions=100L)
kde_gl_gaussian <- kernel_density_estimator(gaussian_kernel, samples, gl_gaussian, subdivisions=100L)
kde_pco_gaussian <- kernel_density_estimator(gaussian_kernel, samples, pco_gaussian, subdivisions=100L)

kde_cv_sigmoid <- kernel_density_estimator(sigmoid_kernel, samples, cv_sigmoid, subdivisions=200L)
kde_gl_sigmoid <- kernel_density_estimator(sigmoid_kernel, samples, gl_sigmoid, subdivisions=200L)
kde_pco_sigmoid <- kernel_density_estimator(sigmoid_kernel, samples, pco_sigmoid, subdivisions=200L)


# comparing the ISE
ISE <- function(kde, dens, subdivisions=100L){
  support <- range(kde$support, dens$support)
  integrate(function(x){(kde$fun(x)-dens$fun(x))^2}, lower=support[1], upper=support[2], subdivisions=subdivisions)$value
}
# ISE for the gaussian kernel
cat("ISE for cross validation on the gaussian kernel: ", ISE(kde_cv_gaussian, dens))
cat("ISE for Goldenshluger-Lepski method on the gaussian kernel: ", ISE(kde_gl_gaussian, dens))
cat("ISE for PCO method on the gaussian kernel: ", ISE(kde_pco_gaussian, dens))

# ISE for the sigmoid kernel
cat("ISE for cross validation method on the sigmoid kernel: ", ISE(kde_cv_sigmoid, dens))
cat("ISE for Goldenshluger-Lepski method method on the sigmoid kernel: ", ISE(kde_gl_sigmoid, dens))
cat("ISE for PCO method method on the sigmoid kernel: ", ISE(kde_pco_sigmoid, dens))

table <- data.frame(method=c("cross_validation", "goldenshluger_lepski", "pco_method"),
                    gaussian=c(ISE(kde_cv_gaussian, dens), ISE(kde_gl_gaussian, dens), ISE(kde_pco_gaussian, dens)),
                    sigmoid=c(ISE(kde_cv_sigmoid, dens), ISE(kde_gl_sigmoid, dens),  ISE(kde_pco_sigmoid, dens)))

# plotting KDE and functionx_lim
x_lim_lower <- -0.5
x_lim_upper<- 1.5
x <- seq(from = x_lim_lower, to = x_lim_upper, length.out=1000)
plot(x, dens$fun(x),
     xlim = c(x_lim_lower, x_lim_upper),
     ylim = c(-0.5, 2.5),
     main = "KDE using the gaussian kernel",
     xlab = "",
     ylab = "",
     col = "dark red",
     type = "l",
     lwd = 2

)
legend("topright", legend = c("density", "samples", "PCO method", "Crossvalidation", "Goldenshluger-Lepski"), col = c("dark red","royal blue", "dark green","violet", "steelblue2"), lty = c(1,1,1,1,1), lwd = c(2,1,2,1,1), cex = 0.75)
lines(x,
      kde_cv_gaussian$fun(x), col = "violet")
lines(x,
      kde_gl_gaussian$fun(x), col = "steelblue2")
lines(x,
      kde_pco_gaussian$fun(x), col = "dark green")
points(samples,
       integer(length(samples)),
       pch = ".",
       col = "blue")


plot(x, dens$fun(x),
     xlim = c(x_lim_lower, x_lim_upper),
     ylim = c(-0.5, 2.5),
     main = "KDE using the sigmoid kernel",
     xlab = "",
     ylab = "",
     col = "dark red",
     type = "l",
     lwd = 2
)
legend("topright", legend = c("density", "samples", "PCO method", "Crossvalidation", "Goldenshluger-Lepski"), col = c("dark red","royal blue", "dark green","violet", "steelblue2"), lty = c(1,1,1,1,1), lwd = c(2,1,2,1,1), cex = 0.75)
lines(x,
      kde_cv_sigmoid$fun(x), col = "violet")
lines(x,
      kde_gl_sigmoid$fun(x), col = "steelblue2")
lines(x,
      kde_pco_sigmoid$fun(x), col = "dark green")
points(samples,
       integer(length(samples)),
       pch = ".",
       col = "blue")
