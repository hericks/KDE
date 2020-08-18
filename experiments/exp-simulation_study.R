#code for simulation study vignette
library(tidyverse)
# Custom density
f_dens <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}
support_density <- c(0,1)
# constructing the Density object
dens <- Density(f_dens, support_density, subdivisions=1000L)

x_lim <- c(dens$support[1] - 0.5, dens$support[2] + 0.5)
grid <- seq(from=x_lim[1], to=x_lim[2], length.out=300)
# plotting the density



plot(grid, dens$fun(grid), xlim = x_lim, ylim=c(-1,2),
     xlab = "",
     ylab = "",
     col = "dark red",
     type = "l",
     lwd = 2)

# setting up a sampler for the density

# Create sampler from custom density
g_den <- Density(dunif, c(0,1))
custom_sampler <- rejection_sampling(dens, g_den, runif, 2)

# Sample from custom sampler
samples <- list()
n_samples <- c(10, 50, 1000)

kernels <- list(rectangular=rectangular, gaussian=gaussian, epanechnikov=epanechnikov)
dens_list <- list(list(dens,custom_sampler), list(Density(dunif,c(0,1)), runif), list(Density(dnorm,c(-15,15)), rnorm))
bandwidth <- 0.1
par(mfrow = c(1, 3))
for(j in seq_along(dens_list)) {
  d <- dens_list[[j]]
    for (i in seq_along(kernels)) {
      for(n in n_samples) {
      samples <- d[[2]](n)
      name <- names(kernels)[[i]]
      k <- kernels[[i]]
        kde <- kernel_density_estimator(k, samples, bandwidth = bandwidth)
        if(j < 3){
          grid <- seq(from=x_lim[1], to=x_lim[2], length.out=300)
          plot(
            grid,
            d[[1]]$fun(grid),
            xlim = x_lim,
            ylim = c(-0.5, 2),
            main = paste(length(samples), "samples"),
            xlab = "",
            ylab = "",
            col = "dark red",
            type = "l",
            lwd = 2
          )
        }
        else{
          plot(
            grid <- seq(from = -16 ,to = 16, length.out=3000),
            d[[1]]$fun(grid),
            xlim = c(-4, 4),
            ylim = c(-0.2, 1),
            main = paste(length(samples), "samples"),
            xlab = "",
            ylab = "",
            col = "dark red",
            type = "l",
            lwd = 2
          )
        }
        lines(grid, kde$fun(grid), col = "orange")
      }
  }
}

# as we can see, as the number of samples increases, the kernel will get more and more irrelevant
# mathematically we will calculate the MISE to get the empirical proof
mise_vec <- c()
reps <- 2

for(j in seq_along(dens_list)) {
  d <- dens_list[[j]]
  for(i in seq_along(kernels)){
    name <- names(kernels)[[i]]
    k <- kernels[[i]]
    for(n in n_samples){
      ise_vec <- c()
      for(rep in 1:reps){
        samples <- d[[2]](n)
        kde <- kernel_density_estimator(k, samples, bandwidth=bandwidth)
        ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
      }
      mise_vec <- c(mise_vec, mean(ise_vec))
    }
  }
}

mise <- array(mise_vec,
              dimnames = list("n_samples" =n_samples, "kernels"=c("rectangular","gaussian", "epanechnikov"),
                              "density"=c("custom density", "uniform distribution", "normal distribution")),
              dim=c(length(n_samples),length(kernels),length(dens_list)))
print(mise)

# -> we will set the kernel to epanechnikov and work on the largest set of samples
# the important parameter for a good estimation is the bandwidth
# above, we did set the bandwidth to h=0.1
# as you can see, the bandwidth is a much more important parameter than the kernel:
par(mfrow = c(1,1))
bandwidth_set <- list(list(0.3, "dark red"), list(0.01, "dark green"), list(0.001, "orange"))
kernel <- epanechnikov
n_samples <- 1000
samples <- custom_sampler(n_samples)
plot(grid, dens$fun(grid), xlim = x_lim, ylim=c(-0.1,2),
     main="comparison of KDE with different bandwidths",
     xlab = "",
     ylab = "",
     type = "l",
     lwd = 2)
legend("topright", title= "bandwidths", legend = c(0.3, 0.01, 0.001), col = c("dark red", "dark green", "orange"), lty = c(1,1,1), lwd = c(1,1,1), cex = 1.2)
for(h in bandwidth_set){
  kde <- kernel_density_estimator(kernel ,samples, bandwidth = h[[1]])
  lines(grid, kde$fun(grid), col = h[[2]])
}

mise_vec <- c()
reps <- 2
for(j in seq_along(dens_list)) {
  d <- dens_list[[j]]
  for(h in bandwidth_set){
    ise_vec <- c()
    for(rep in 1:reps){
      samples <- d[[2]](n_samples)
      kde <- kernel_density_estimator(kernel, samples, bandwidth=h[[1]])
      ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
    }
    mise_vec <- c(mise_vec, mean(ise_vec))
  }
}

mise <- matrix(mise_vec, nrow=c(length(bandwidth_set)),
               dimnames = list("bandwidhts"= c(0.3, 0.01, 0.001),
                               "density"=c("custom density", "uniform distribution", "normal distribution")))
print(mise)


# TODO: 1. bandweitenschätzer vorstellen (mathematisch?),
# sweeter plot f?r defaultwerte auf unserer custom_density:
# TODO: main und legend mit reinwerfen (als listen?)
plot_comparison(show_diff=FALSE, reps=2)

# 2. lambda/kappa wählen
# hence we estimate a upper bound for the variance, we have tuning parameters for the pco_method and goldenshluger_lepski

# first, lets take a look at goldenshluger_lepski and the kappa parameter
# in literature, they set kappa=1.2

# kappa_set <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
ns <- 100
reps <- 2
kappa_set <- list(1, 2)
dens_list <- list(custom_dens=dens, dunif=Density(dunif,c(0,1)), dnorm=Density(dnorm,c(-15,15)))
sampler_list <- list(custom_sampler, runif, rnorm)
kernel_list <- list(epanechnikov=epanechnikov)
ise_kappa <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, kernels=kernel_list, kappa_set=kappa_set, ns = ns, reps = reps)
mise_kappa <- calculate_mise(ise_kappa)

# bestes kappa wird hier ... sein

# now we try to tune our lambda
# lambda_set <- c(1, 1.2, 1.4, 1.6, 1.8, 2)
ns <- 100
reps <- 2
lambda_set <- list(1,2)
dens_list <- list(custom_dens=dens, dunif=Density(dunif,c(0,1)), dnorm=Density(dnorm,c(-15,15)))
sampler_list <- list(custom_sampler, runif, rnorm)
kernel_list <- list(epanechnikov=epanechnikov)
ise_lambda <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, kernels=kernel_list, lambda_set=lambda_set, ns = ns, reps = reps)
mise_lambda <- calculate_mise(ise_lambda)

# bestes lambda wird hier ... sein

# 3. Schätzer mit optimalen lambda/kappa
# nun vergleichen wir die bandweitensch?tzer untereinander
# lambda, kappa fest wie oben gew?hlt

kappa_set <- list(1.2)
lambda_set <- list(1)
reps <- 2
ns <- c(1000)

# mehr densities?
# als erstes schauen wir, welcher Bandweitensch?tzer bei 1000 samples auf unseren 3 densities am besten performen w?rde
for(i in seq_along(dens_list)){
  d <- dens_list[[i]]
  xlim <-  c(d$support[1] - 1, d$support[2] + 1)
  plot_comparison(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim[1], xlim_upper=xlim[2], reps=reps, ns=ns)
}

ise <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, reps=reps,ns=ns)
mise_ns_comp <- calculate_mise(ise)
mise_ns_comp %>%
  group_by(n, bandwidth_estimators) %>%
  summarize(mean_mise=mean(mise), mean_sd_ise=mean(sd_ise))


# TODO: performance
# we will make a small performance comparison

# da nicht immer sehr viele daten vorhanden sind, betrachten wir nun, welcher Bandweitensch?tzer sich bei verschiedenen sample sizes am besten verh?lt
c(10, 50, 100, 1000)
par(mfrow=c(1,4))
for(i in seq_along(dens_list)){
  par(mfrow=c(1,4))
  d <- dens_list[[i]]
  xlim <-  c(d$support[1] - 1, d$support[2] + 1)
  plot_comparison(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim[1], xlim_upper=xlim[2], reps=reps, ns=ns)
}

par(mfrow=c(1,1))

ise <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, reps=reps, ns=c(10, 50, 100, 1000))
mise_ns_comp <- calculate_mise(ise)
mise_ns_comp %>%
  group_by(n, bandwidth_estimators) %>%
  summarize(mean_mise=mean(mise), mean_sd_ise=mean(sd_ise))


# 4. MSE am Rand von Dichte mit kompaktem support?

#microbenchmark(expr1, expr2, times=25L)
# TODO: 1. performance vergleich, 1.5 legenden etc 2. wie kann man objekte abspeichern und in vignette einbinden? 3. welche parameter bei welchem test?

