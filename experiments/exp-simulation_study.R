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

#TODO: Iteration fixen, das Problem ist dass bei #Problem immer d das dritte Element in dens_list ist aber in der 2.
#Schleife wo die samples erstellt werden, läuft die iteration über dens_list richtig (siehe Plot mit KDE)

#hinzugefügt: iteration über Densities! hier ist noch was nciht ganz fertig!
kernels <- list(rectangular=rectangular, gaussian=gaussian, epanechnikov=epanechnikov)
dens_list <- list(list(dens,custom_sampler), list(Density(dunif,c(0,1)), runif), list(Density(dnorm,c(-15,15)), rnorm))
for(d in dens_list) {
  for (n in n_samples) {
    samples <- c(samples, list(d[[2]](n)))   #das funktioniert
    for (i in seq_along(kernels)) {
      name <- names(kernels)[[i]]
      k <- kernels[[i]]
      par(mfrow = c(1, 3))
      for (s in samples) {
        kde <- kernel_density_estimator(k, s, bandwidth = 0.1)
        if(identical(d[[2]], rnorm)) {                             #Problem: hier wird nicht mehr über d iteriert..
          grid <- seq(from = -16 ,to = 16, length.out=3000)
          plot(
            grid <- seq(from = -16 ,to = 16, length.out=3000),
            d[[1]]$fun(grid),
            xlim = c(-4, 4),
            ylim = c(-0.2, 1),
            main = paste(length(s), "samples"),
            xlab = "",
            ylab = "",
            col = "dark red",
            type = "l",
            lwd = 2
          )
        }
        else{
        plot(
          grid,
          d[[1]]$fun(grid),
          xlim = x_lim,
          ylim = c(-0.5, 2),
          main = paste(length(s), "samples"),
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
}

# as we can see, as the number of samples increases, the kernel will get more and more irrelevant
# mathematically we will calculate the MISE to get the empirical proof
mise_vec <- c()
reps <- 100
# TODO: hier noch die anderen denseties einfügen für den mise

for(i in seq_along(kernels)){
  name <- names(kernels)[[i]]
  k <- kernels[[i]]
  for(n in n_samples){
    print(i)
    ise_vec <- c()
    for(rep in 1:reps){
      samples <- custom_sampler(n)
      kde <- kernel_density_estimator(k, samples, bandwidth=0.1)
      ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
    }
    mise_vec <- c(mise_vec, mean(ise_vec))
  }
}

mise <- matrix(mise_vec, nrow=length(n_samples), dimnames = list("n_samples" =c(10,50,1000), "kernels"=c("rectangular","gaussian", "epanechnikov")))
print(mise)

# -> we will set the kernel to gaussian and work on the largest set of samples
# the important parameter for a good estimation is the bandwidth
# above, we did set the bandwidth to h=0.1
# as you can see, the bandwidth is a much more important parameter than the kernel:
par(mfrow = c(1,1))
bandwidth_set <- list(list(0.5,"dark green"), list(0.1, "dark red"), list(0.01, "steelblue3"), list(0.001, "orange"))
kernel <- gaussian
n_samples <- 1000
samples <- custom_sampler(n_samples)
plot(grid, dens$fun(grid), xlim = x_lim, ylim=c(-1,2),
     main="comparison of KDE with different bandwidths",
     xlab = "",
     ylab = "",
     type = "l",
     lwd = 2)
legend("topright", title= "bandwidths", legend = c(0.5, 0.1, 0.01, 0.001), col = c("dark green","dark red", "steelblue3", "orange"), lty = c(1,1,1), lwd = c(1,1,1), cex = 1.2)
for(h in bandwidth_set){
  kde <- kernel_density_estimator(kernel ,samples, bandwidth = h[[1]])
  lines(grid, kde$fun(grid), col = h[[2]])
}

mise_vec <- c()
reps <- 10
for(h in bandwidth_set){
  ise_vec <- c()
  for(rep in 1:reps){
    samples <- custom_sampler(n_samples)
    kde <- kernel_density_estimator(kernel, samples, bandwidth=h[[1]])
    ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
  }
  mise_vec <- c(mise_vec, mean(ise_vec))
}
mise <- matrix(mise_vec, nrow=length(bandwidth_set), dimnames = list("bandwidhts"= c(0.5, 0.1, 0.01, 0.001)))
print(mise)


# TODO: 1. bandweitenschätzer vorstellen (mathematisch?),
# 2. lambda/kappa wählen
# 3. Schätzer mit optimalen lambda/kappa
