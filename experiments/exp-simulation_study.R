dens_list <- list(
  uniform_distri=Density(dunif, c(0,1)),
  normal_distri=Density(dnorm, c(-15,15))

)

sampler_list <- list(
  uniform_distri=runif,
  normal_distri=rnorm
)

#data_ise <- compare_ise(dens_list, sampler_list)
#mse_data <- calculate_mise(data_ise)
#mse_data

#code for simulation study vignette

# Custom density
f_dens <- function(x) {
  ret <- 1 + sin(2*pi*x)
  ret[x < 0 | 1 < x] <- 0
  ret
}
support_density <- c(0,1)
# constructing the Density object
dens <- Density(f_dens, support_density, subdivisions=1000L)

x_lim <- c(dens$support[1] - 1, dens$support[2] + 1)
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
for(n in n_samples){
  samples <- c(samples, list(custom_sampler(n)))
}

kernels <- list(rectangular=rectangular, gaussian=gaussian, epanechnikov=epanechnikov)

for(i in seq_along(kernels)){
  name <- names(kernels)[[i]]
  k <- kernels[[i]]
  par(mfrow=c(1,3))
  for(s in samples){
  kde <- kernel_density_estimator(k,s, bandwidth=0.1)
  plot(grid, dens$fun(grid), xlim = x_lim, ylim=c(-1,2),
       main=paste(length(s), "samples"),
       xlab = "",
       ylab = "",
       col = "dark red",
       type = "l",
       lwd = 2)
  lines(grid, kde$fun(grid))
  }
}

# as we can see, as the number of samples increases, the kernel will get more and more irrelevant
# mathematically we will calculate the MISE to get the empirical proof
mise_vec <- c()
reps <- 100
# TODO: bandweite, die das gut demonstriert wählen (auch in plots)
for(i in seq_along(kernels)){
  name <- names(kernels)[[i]]
  k <- kernels[[i]]
  for(n in n_samples){
    print(i)
    ise_vec <- c()
    for(rep in 1:reps){
      samples <- custom_sampler(n)
      kde <- kernel_density_estimator(k, samples, bandwidth=0.5)
      ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
    }
    mise_vec <- c(mise_vec, mean(ise_vec))
  }
}

#TODO: benennen der spalten und reihen
mise <- matrix(mise_vec, nrow=length(n_samples))
mise

# -> we will set the kernel to gaussian and work on the largest set of samples
# the important parameter for a good estimation is the bandwidth
# above, we did set the bandwidth to h=0.1
# as you can see, the bandwidth is a much more important parameter than the kernel:
# TODO: legende und verschiedene Farben für bandweiten
par(mfrow = c(1,1))
bandwidth_set <- list(1, 0.5, 0.1, 0.01, 0.001)
kernel <- gaussian
n_samples <- 1000
samples <- custom_sampler(n_samples)
plot(grid, dens$fun(grid), xlim = x_lim, ylim=c(-1,2),
     main="bandwidths",
     xlab = "",
     ylab = "",
     col = "dark red",
     type = "l",
     lwd = 2)
for(h in bandwidth_set){
  kde <- kernel_density_estimator(kernel ,samples, bandwidth=h)
  lines(grid, kde$fun(grid))
}

mise_vec <- c()
reps <- 10
for(h in bandwidth_set){
  ise_vec <- c()
  for(rep in 1:reps){
    samples <- custom_sampler(n_samples)
    kde <- kernel_density_estimator(kernel, samples, bandwidth=h)
    ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
  }
  mise_vec <- c(mise_vec, mean(ise_vec))
}
#TODO: benennen der spalten und reihen
mise <- matrix(mise_vec, ncol=length(bandwidth_set))
print(mise)


# TODO: 1. bandweitenschätzer vorstellen (mathematisch?),
# 2. lambda/kappa wählen
# 3. Schätzer mit optimalen lambda/kappa
