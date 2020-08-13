# print integrable function?????????????

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

# using a builtin Kernel
gaussian_kernel <- gaussian
print(gaussian_kernel)

# or build a Kernel object yourself: this is the epanechnikov function
f_ker <- function(x){
  return(3/4 * (1 - x^2) * (abs(x) <= 1))
}
support_epanechnikov <- c(-1,1)
# constructing the Kernel object
epanechnikov_kernel <- Kernel(f_ker, support_epanechnikov, subdivisions=10L)
print(epanechnikov_kernel)


# Create sampler from custom density
g_den <- Density(dunif, c(0,1))
custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)
# Sample from custom sampler
num_samples <- 25
samples <- custom_sampler(num_samples)

# setting up bandwidth sets
bandwidth_set <- logarithmic_bandwidth_set(from=1/length(samples), to=1, length.out=10)

# bandwidth estimation
# cross-validation method
cv_gaussian <- cross_validation(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
cv_epanechnikov <- cross_validation(epanechnikov_kernel, samples, bandwidths=bandwidth_set, subdivisions=225L)
cat("bandwidth for gaussian kernel: ", cv_gaussian, "\nbandwidth for epanechnikov kernel: ", cv_epanechnikov)

#  goldenshluger-lepski method
gl_gaussian <- goldenshluger_lepski(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
gl_epanechnikov <- goldenshluger_lepski(epanechnikov_kernel, samples, bandwidths=bandwidth_set, subdivisions=200L)
cat("bandwidth for gaussian kernel: ", gl_gaussian, "\nbandwidth for epanechnikov kernel: ", gl_epanechnikov)

#  Penalized Comparison to Overfitting
pco_gaussian <- pco_method(gaussian_kernel, samples, bandwidths=bandwidth_set, subdivisions=100L)
pco_epanechnikov <- pco_method(epanechnikov_kernel, samples, bandwidths=bandwidth_set, subdivisions=200L)
cat("bandwidth for gaussian kernel: ", pco_gaussian, "\nbandwidth for epanechnikov kernel: ", pco_epanechnikov)


# Kernel Density Estimation


kde_cv_gaussian <- kernel_density_estimator(gaussian_kernel, samples, cv_gaussian, subdivisions=100L)
kde_gl_gaussian <- kernel_density_estimator(gaussian_kernel, samples, gl_gaussian, subdivisions=100L)
kde_pco_gaussian <- kernel_density_estimator(gaussian_kernel, samples, pco_gaussian, subdivisions=100L)



# plotting KDE and functionx_lim
x_lim_lower <- -0.5
x_lim_upper<- 1.5
x <- seq(from = x_lim_lower, to = x_lim_upper, length.out=1000)
plot(x, dens$fun(x),
     xlim = c(x_lim_lower, x_lim_upper),
     ylim = c(-0.5, 2),
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
     ylim = c(-0.5, 2),
     main = "KDE using the epanechnikov kernel",
     xlab = "",
     ylab = "",
     col = "dark red",
     type = "l",
     lwd = 2
)
legend("topright", legend = c("density", "samples", "PCO method", "Crossvalidation", "Goldenshluger-Lepski"), col = c("dark red","royal blue", "dark green","violet", "steelblue2"), lty = c(1,1,1,1,1), lwd = c(2,1,2,1,1), cex = 0.75)
lines(x,
      kde_cv_epanechnikov$fun(x), col = "violet")
lines(x,
      kde_gl_epanechnikov$fun(x), col = "steelblue2")
lines(x,
      kde_pco_epanechnikov$fun(x), col = "dark green")
points(samples,
       integer(length(samples)),
       pch = ".",
       col = "blue")
