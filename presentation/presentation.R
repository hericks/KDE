# Setup
dev.off()
quartz(title="Plots",
       width = 5.25,
       height= 3.75,
       dpi=120)

par(mar=double(4))
plot(1, type="n", xlab="", ylab="", yaxt="n", xaxt="n")
options(viewer=NULL)

### Introduction to kernel density estimation
set.seed(17)
samples <- sort(rnorm(10))

# empty plot
par(mar=c(2.5, 0, 0, 0))
clear_plot <- function()
  plot(1, type="n", xlab="", ylab="", xlim=c(-3, 3), ylim=c(0, 0.7))
clear_plot()

# plot samples
plot_samples <- function()
  points(samples, double(length(samples)) - 0.015, col="blue", pch=".", cex=2)
plot_samples()

# plot density
grid <- seq(-3, 3, length.out=1000L)
plot_density <- function() lines(grid, dnorm(grid), col="red", lty="dashed")
plot_density()

# kernel
kernel <- function(grid) (abs(grid) <= 1)/2

# plot estimator given first sample/second sample #
lines(grid, kernel(grid - samples[1]))

lines(grid, kernel(grid - samples[2]))

## plot mean
clear_plot()
plot_density()
lines(grid, (kernel(grid - samples[1]) + kernel(grid - samples[2]))/2)
plot_samples()

##
clear_plot()
plot_density()
plot_samples()
for (sample in samples) {
  lines(grid, kernel(grid - sample)/length(samples))
}

#
estimator_up_to <- function(n) {
  function(x) {
    res <- 0
    for (i in 1:n) {
      res <- res + kernel(x - samples[i])
    }
    res/length(samples)
  }
}

lines(grid, estimator_up_to(2)(grid))
lines(grid, estimator_up_to(3)(grid))
lines(grid, estimator_up_to(4)(grid))
lines(grid, estimator_up_to(10)(grid))

## plot without isolated kernels
clear_plot()
plot_density()
lines(grid, estimator_up_to(10)(grid))
plot_samples()

# example for different kernels #
triangular <- function(x) {
  return((1 - abs(x)) * (abs(x) <= 1))
}

epanechnikov <- function(x){
  return(3/4 * (1 - x^2) * (abs(x) <= 1))
}

gaussian <- function(x){
  return(1/sqrt(2 * pi) * exp(-1/2 * x^2))
}

grid <- seq(-2.5, 2.5, length.out=300)

par(mfrow=c(2,2), mar=c(0, 0, 0, 0))
plot(grid, kernel(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, triangular(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, epanechnikov(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, gaussian(grid), type="l",  yaxt = "n", xaxt = "n")
rm(triangular, epanechnikov, gaussian)


### Introduction to the KDE package #
library(KDE)

### Object Structure

# integrable functions

# numeric compact support
IntegrableFunction(dnorm, support=c(-Inf, Inf))
all.equal(dnorm(-15), 0)

f <- IntegrableFunction(dnorm, support=c(-15, 15))
f <- IntegrableFunction(dnorm, mean=2)
str(f)

# subdivisions parameter
integrate_primitive(dnorm, lower=-15, upper=15, subdivisions = 25L)$value
integrate_primitive(dnorm, lower=-15, upper=15, subdivisions = 1000L)$value

f
rm(f)


### Kernels and densities #

# kernels
Kernel(dnorm, subdivisions = 25L)
Kernel(dnorm, subdivisions = 1000L)

Kernel(function(x) (abs(x) <= 1)/2, support=c(-1, 1))
rectangular

## plot built-in kernels
par(mfrow=c(2,2), mar=c(0, 0, 0, 0))
grid <- seq(-2.5, 2.5, length.out=1000)
plot(grid, rectangular$fun(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, triangular$fun(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, epanechnikov$fun(grid), type="l",  yaxt = "n", xaxt = "n")
plot(grid, gaussian$fun(grid), type="l",  yaxt = "n", xaxt = "n")

# densities #
Density(dnorm)



### Construction of kernel density estimator
set.seed(17)

samples <- rnorm(10)

## Recreate setting
grid <- seq(-3, 3, length.out=300L)
par(mfrow=c(1,1), mar=c(2.5, 0, 0, 0))
rebuild_plot <- function() {
  par(mar=c(2.5, 0, 0, 0))
  plot(1, type="n", xlab="", ylab="", xlim=c(-3, 3), ylim=c(0, 0.7))
  points(samples, double(length(samples)) - 0.015, col="blue", pch=".", cex=3)
  lines(grid, dnorm(grid), col="red", lty="dashed")
}
rebuild_plot()


kde <- kernel_density_estimator(rectangular, samples, 1)
lines(grid, kde$fun(grid), type="l")

rebuild_plot()
kde <- kernel_density_estimator(cosine, samples, 1)
lines(grid, kde$fun(grid), type="l")

# bandwidth comparison
rebuild_plot()
kde_underfit <- kernel_density_estimator(gaussian, samples, 1)
lines(grid, kde_underfit$fun(grid), type="l", col="grey")
kde_overfit <- kernel_density_estimator(gaussian, samples, 0.15)
lines(grid, kde_overfit$fun(grid), type="l", col="grey")
kde <- kernel_density_estimator(gaussian, samples, 0.6)
lines(grid, kde$fun(grid), type="l")

estimators <- list(large_bw=kde_underfit, medium_bw=kde, small_bw=kde_overfit)
ises <- sapply(estimators, function(estimator) {
  integrate_primitive(function(x) (estimator$fun(grid) - dnorm(grid))^2,
                      estimator$support[1],
                      estimator$support[2],
                      subdivisions = 5000L)$value
})
ises



### Bandwidth estimation #
# open vignette

set.seed(17)
samples <- rnorm(500)

##
clear_plot()
plot_samples()
plot_density()

kernel <- gaussian
bandwidths <- logarithmic_bandwidth_set(1/500, 1, 15)

cv_bandwidth <- cross_validation(kernel, samples, bandwidths)
cv_bandwidth
cv_estimator <- kernel_density_estimator(kernel, samples, cv_bandwidth)
lines(grid, cv_estimator$fun(grid), col="darkorchid2")

gl_bandwidth <- goldenshluger_lepski(kernel, samples, bandwidths)
gl_bandwidth
gl_estimator <- kernel_density_estimator(kernel, samples, gl_bandwidth)
lines(grid, gl_estimator$fun(grid), col="chartreuse")

pco_bandwidth <- pco_method(kernel, samples, bandwidths)
pco_bandwidth
pco_estimator <- kernel_density_estimator(kernel, samples, pco_bandwidth)
lines(grid, pco_estimator$fun(grid), col="deepskyblue")

legend("topright",
       legend=c("density", "cross-validation", "goldenshluger-lepski", "pco"),
       lty = c(2, 1, 1, 1),
       col = c("red", "darkorchid2", "chartreuse", "deepskyblue"),
       text.width = 2,
       cex = 0.75)

estimators <- list(cv=cv_estimator, gl=gl_estimator, pco=pco_estimator)
ises <- sapply(estimators, function(estimator) {
  integrate_primitive(function(x) (estimator$fun(grid) - dnorm(grid))^2,
                      estimator$support[1],
                      estimator$support[2],
                      subdivisions = 5000L)$value
})
ises



### Shiny
shiny_kde()



### Obstacles: numeric integration
shifted_kernel <- kernel_transform(rectangular, 10, 0.01)

grid <- seq(0, 20, by=0.01)
plot(grid, shifted_kernel$fun(grid), type="l")

integrate(shifted_kernel$fun, -Inf, Inf)
integrate(shifted_kernel$fun, shifted_kernel$support[1], shifted_kernel$support[2])

estimator <- kernel_density_estimator(rectangular, c(-10, -5, 10), 0.01)

grid <- seq(-12, 12, by=0.0001)
plot(grid, estimator$fun(grid), type="l")

integrate(estimator$fun, estimator$support[1], estimator$support[2])
estimator$support

integrate_primitive(estimator$fun, estimator$support[1], estimator$support[2])$value


### Outlook

# convergence
integrate_primitive(dnorm, -10, 10, subdivisions = 500, check = TRUE)
integrate_primitive(dnorm, -10, 10, subdivisions = 1000, check = TRUE)
integrate_primitive(dnorm, -10, 10, subdivisions = 5000, check = TRUE)
integrate_primitive(dnorm, -10, 10, subdivisions = 10000, check = TRUE)

# divergence
integrate_primitive(function(x) 1/x, 0, 1, subdivisions = 500, check = TRUE)
integrate_primitive(function(x) 1/x, 0, 1, subdivisions = 1000, check = TRUE)
integrate_primitive(function(x) 1/x, 0, 1, subdivisions = 5000, check = TRUE)
integrate_primitive(function(x) 1/x, 0, 1, subdivisions = 10000, check = TRUE)
