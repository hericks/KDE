library(tidyverse)

# TODO::length.out hochsetzen, length.out = 10
compare <- function(eval_points,
                    funs = list(runif),
                    bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                    ns = 1000,
                    kernels = list(gaussian),
                    lambda_set = list(1),
                    kappa_set = list(1.2),
                    reps = 400,
                    length.out = 10) {
  if (!is.list(kernels))
    kernels <- as.list(kernels)
  if (!is.list(ns))
    ns <- as.list(ns)

  if(is.list(eval_points)){
    if(!(length(unique(sapply(eval_points, length))) == 1)) stop("same length for all eval_points vectors are needed")
    if(length(eval_points) != length(funs) && length(eval_points) > 1) stop("eval_point set has to be of same length as funs")
  }

  if (length(kappa_set) > 1 &&
      length(lambda_set) > 1) {
    stop("how is this supposed to look like?")
  }
  else if (length(kappa_set) > 1) {
    bandwidth_estimators <- list(goldenshluger_lepski=goldenshluger_lepski)  }
  else if (length(lambda_set) > 1) {
    bandwidth_estimators <- list(pco_method=pco_method)  }


  m <- length(funs) * length(ns) * length(kernels)
  m <-
    m * (any(sapply(
      bandwidth_estimators, identical, cross_validation
    ))
    + length(kappa_set) * any(
      sapply(bandwidth_estimators, identical, goldenshluger_lepski)
    )
    + length(lambda_set) * any(sapply(
      bandwidth_estimators, identical, pco_method
    )))
  res <- array(NA, dim = c(length(eval_points[[1]]), reps, m))
  cnt <- 1
  subdivisions <- 1000L
  for (i in seq_along(funs)) {
    if(length(eval_points) > 1){x_grid <- eval_points[[i]]
    } else{
      if(is.list(eval_points)) {x_grid <- eval_points[[1]]
      }else{
        x_grid <- eval_points
      }
    }
    f <- funs[[i]]
    for (n in ns) {
      bandwidths <-
        logarithmic_bandwidth_set(from = 1 / n,
                                  to = 1,
                                  length.out = length.out)
      for (k in kernels) {
        for (est in bandwidth_estimators) {
          cat("fun number:", i)
          if (identical(est, goldenshluger_lepski)) {
            for (kappa in kappa_set) {
              res[, , cnt] <- replicate(reps, {
                samples <- f(n)
                bandwidth <-
                  goldenshluger_lepski(k, samples, bandwidths, kappa, subdivisions)

                kde <-
                  kernel_density_estimator(k, samples, bandwidth, subdivisions)
                if (any(kde$fun(x_grid) < 0)) {
                  print(kde)
                }
                kde$fun(x_grid)
              })
              cnt <- cnt + 1
            }
          } else if (identical(est, pco_method)) {
            for (lambda in lambda_set) {
              res[, , cnt] <- replicate(reps, {
                samples <- f(n)

                bandwidth <-
                  pco_method(k, samples, bandwidths, lambda, subdivisions)

                kde <-
                  kernel_density_estimator(k, samples, bandwidth, subdivisions)
                if (any(kde$fun(x_grid) < 0)) {
                  print(kde)
                }

                kde$fun(x_grid)
              })
              cnt <- cnt + 1
            }
          } else if (identical(est, cross_validation)) {
            res[, , cnt] <- replicate(reps, {
              samples <- f(n)

              bandwidth <-
                cross_validation(k, samples, bandwidths, subdivisions)

              kde <-
                kernel_density_estimator(k, samples, bandwidth, subdivisions)
              if (any(kde$fun(x_grid) < 0)) {
                print(kde)
              }

              kde$fun(x_grid)
            })
            cnt <- cnt + 1
          }
        }
      }

    }
  }
  res
}


plot_with_confidence_band <- function(x, Y, col) {
  rgb <- col2rgb(col) / 255
  col_alpha <- rgb(rgb[1], rgb[2], rgb[3], 0.2)
  v <- apply(Y, 1, function(x)
    c(mean(x), sd(x)))
  lines(x, v[1, ], lwd = 2, col = col)
  polygon(c(x, rev(x)),
          c(v[1, ] + v[2, ], rev(v[1,] - v[2,])),
          col = col_alpha,
          border = NA)
  lines(x, v[1,] + v[2,], lwd = 1, col = col)
  lines(x, v[1,] - v[2,], lwd = 1, col = col)
}

plot_comparison_objects <- function(dens = Density(dunif, c(0, 1)),
                                    dens_sampler = runif,
                                    xlim_lower = -1,
                                    xlim_upper = 2,
                                    main = NA,
                                    legend = NULL,
                                    show_diff = TRUE,
                                    split = TRUE,
                                    bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                                    reps = 2,
                                    kappa = list(1.2),
                                    lambda = list(1),
                                    ...) {
  dens_fun <- dens$fun
  x_grid <- seq(xlim_lower, xlim_upper, length.out = 300)

  dens_eval <- dens_fun(x_grid)
  eval_points <- list(x_grid)
  if (length(kappa) > 1 &&
      length(lambda) > 1) {
    stop("how is this plot supposed to look like?")
  }
  else if (length(kappa) > 1) {
    res <-
      compare(
        eval_points,
        list(dens_sampler),
        reps = reps,
        bandwidth_estimators = list(goldenshluger_lepski),
        kappa_set = kappa,
        ...
      )
  }
  else if (length(lambda) > 1) {
    res <-
      compare(
        eval_points,
        list(dens_sampler),
        reps = reps,
        bandwidth_estimators = list(pco_method),
        lambda_set = lambda,
        ...
      )
  } else {
    res <-
      compare(
        eval_points,
        list(dens_sampler),
        reps = reps,
        bandwidth_estimators = bandwidth_estimators,
        ...
      )
  }
  m <- dim(res)[3]
  bw_len <- length(bandwidth_estimators)
  if (isTRUE(split) && bw_len > 1 && m > 3) {
    m_bw_len <- m / bw_len
    ret <- list()
    for (i in 1:m_bw_len) {
      res_sub <- res[, , (1 + (i - 1) * bw_len):(i * bw_len)]
      ret <- c(ret, list(list(res_sub, dens_fun, x_grid, dens_eval)))
    }
  }else{
    ret <- list(res, dens_fun, x_grid, dens_eval)
  }
  ret
}


plot_comparison <- function(dens = Density(dunif, c(0, 1)),
                            dens_sampler = runif,
                            xlim_lower = -1,
                            xlim_upper = 2,
                            ylim_lower=NULL,
                            ylim_upper=NULL,
                            main = NA,
                            legend = NULL,
                            show_diff = TRUE,
                            split = TRUE,
                            bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                            reps = 4,
                            kappa = list(1.2),
                            lambda = list(1),
                            objects=NULL,
                            ...) {


  if(!is.null(objects)){
    res <- objects[[1]]
    dens_fun <- objects[[2]]
    x_grid <- objects[[3]]
    dens_eval <- objects[[4]]
  }else{
    objects <- plot_comparison_objects(dens,
                                       dens_sampler,
                                       xlim_lower,
                                       xlim_upper,
                                       main,
                                       legend,
                                       show_diff,
                                       split,
                                       bandwidth_estimators,
                                       reps,
                                       kappa,
                                       lambda,
                                       ...)
    res <- objects[[1]]
    dens_fun <- objects[[2]]
    x_grid <- objects[[3]]
    dens_eval <- objects[[4]]
  }

  m <- dim(res)[3]
  bw_len <- length(bandwidth_estimators)

  if (isTRUE(split) && bw_len > 1 && m > 3) {
    m_bw_len <- m / bw_len
    for (i in 1:m_bw_len) {
      res_sub <- res[, , (1 + (i - 1) * bw_len):(i * bw_len)]
      del <- max(apply(res_sub, c(1, 3), sd))

      # plot results
      par(
        mar = c(0, 0, 1, 0),
        ann = FALSE,
        xaxt = "n",
        yaxt = "n"
      )
      if(is.null(ylim_lower)){
        ylim_lower <- (range(dens_eval) + del * c(-1, 1))[1]
      }
      if(is.null(ylim_upper)){
        ylim_upper <- (range(dens_eval) + del * c(-1, 1))[2]
      }
      ylim <- c(ylim_lower,ylim_upper)
      plot(
        x_grid,
        dens_eval,
        type = "l",
        lwd = 2,
        col = 1,
        ylim = ylim
      )
      grid()
      for (i in 1:bw_len) {
        plot_with_confidence_band(x_grid, res_sub[, , i], col = i + 1)
      }
      title(main = main)
      if (!is.null(legend))
        legend(
          "topright",
          legend = legend,
          lwd = 2,
          col = 2:(bw_len + 1)
        )

      # second plot (difference between true f and estimation)
      if (show_diff) {
        diff <- res_sub - dens_eval
        plot(
          c(0, 1),
          c(0, 0),
          type = "l",
          lwd = 2,
          col = 1,
          ylim = c(-del, del)
        )
        grid()
        for (i in 1:bw_len)
          plot_with_confidence_band(x_grid, diff[, , i], col = i + 1)
      }
    }
  } else{
    m <- dim(res)[3]
    del <- max(apply(res, c(1, 3), sd))
    # plot results
    par(
      mar = c(0, 0, 1, 0),
      ann = FALSE,
      xaxt = "n",
      yaxt = "n"
    )
    if(is.null(ylim_lower)){
      ylim_lower <- (range(dens_eval) + del * c(-1, 1))[1]
    }
    if(is.null(ylim_upper)){
      ylim_upper <- (range(dens_eval) + del * c(-1, 1))[2]
    }
    ylim <- c(ylim_lower, ylim_upper)
    plot(
      x_grid,
      dens_eval,
      type = "l",
      lwd = 2,
      col = 1,
      ylim = ylim
    )
    grid()
    for (i in 1:m) {
      plot_with_confidence_band(x_grid, res[, , i], col = i + 1)
    }
    title(main = main)
    if (!is.null(legend))
      legend("topright",
             legend = legend,
             lwd = 2,
             col = 2:(m + 1))

    # second plot (difference between true f and estimation)
    if (show_diff) {
      diff <- res - dens_eval
      plot(
        c(0, 1),
        c(0, 0),
        type = "l",
        lwd = 2,
        col = 1,
        ylim = c(-del, del)
      )
      grid()
      for (i in 1:m)
        plot_with_confidence_band(x_grid, diff[, , i], col = i + 1)
    }
  }
}

compare_ise <- function(dens_list = list(dunif=Density(dunif, c(0, 1))),
                        dens_sampler_list = list(runif=runif),
                        bandwidth_estimators = list(cross_validation=cross_validation, goldenshluger_lepski=goldenshluger_lepski, pco_method=pco_method),
                        ns = 50,
                        kernels = list(gaussian=gaussian),
                        lambda_set = list(1),
                        kappa_set = list(1.2),
                        reps = 3,
                        num_eval_points= 300,
                        n_bandwidths = 10){

  if (length(kappa_set) > 1 &&
      length(lambda_set) > 1) {
    stop("how is this supposed to look like?")
  }
  else if (length(kappa_set) > 1) {
    bandwidth_estimators <- list(goldenshluger_lepski=goldenshluger_lepski)  }
  else if (length(lambda_set) > 1) {
    bandwidth_estimators <- list(pco_method=pco_method)}

  eval_points_set <- list()
  for(d in dens_list){
    eval_points_range <-  c(d$support[1] - 1, d$support[2] + 1)
    eval_points <- seq(from=eval_points_range[1],to=eval_points_range[2], length.out=num_eval_points)
    eval_points <- list(eval_points)
    eval_points_set <- c(eval_points_set, eval_points)
  }
  time <-system.time(res <- compare(eval_points=eval_points_set, funs=dens_sampler_list, ns=ns, kernels=kernels,
                                   bandwidth_estimators=bandwidth_estimators,
                                    lambda_set=lambda_set, kappa_set=kappa_set, reps=reps, length.out=n_bandwidths))
  print(time)

  f_true <- array(NA, dim=c(length(eval_points_set[[1]]), length(dens_list)))
  for(i in seq_along(dens_list)){
    f <- dens_list[[i]]$fun
    f_true[,i] <- f(eval_points_set[[i]])
  }
  m <- length(ns) * length(kernels)
  m <-
    m * (any(sapply(
      bandwidth_estimators, identical, cross_validation
    ))
    + length(kappa_set) * any(
      sapply(bandwidth_estimators, identical, goldenshluger_lepski)
    )
    + length(lambda_set) * any(sapply(
      bandwidth_estimators, identical, pco_method
    )))

  f_true <- f_true[,rep(seq_along(dens_list), each=m)]
  diff <-array(NA, dim=dim(res))
  for(j in 1:dim(res)[3]) diff[,,j] <- res[,,j] - f_true[,j]
  diff_sq <-apply(diff^2,c(2, 3), sum)
  ise <- (diff_sq * (eval_points_range[2] - eval_points_range[1])) / num_eval_points

  if (length(kappa_set) > 1 &&
      length(lambda_set) > 1) {
    stop("how is this supposed to look like?")
  }
  else if (length(kappa_set) > 1) {
    opts <- as_tibble(expand.grid(kappa=kappa_set, bandwidth_estimators=names(bandwidth_estimators), kernel=names(kernels), n=ns, den=names(dens_list)))
  }
  else if (length(lambda_set) > 1) {
    opts <-as_tibble(expand.grid(lambda=lambda_set, bandwidth_estimators=names(bandwidth_estimators), kernel=names(kernels), n=ns, den=names(dens_list)))

  } else {
    opts <-as_tibble(expand.grid(bandwidth_estimators=names(bandwidth_estimators), kernel=names(kernels), n=ns, den=names(dens_list)))
  }
  add_column(opts[rep(1:nrow(opts), each=reps), ], ise=as.vector(ise))
}

calculate_mise <- function(data_ise){
  if("kappa" %in% names(data_ise)){
    data_ise %>%
      group_by(den, n, kernel, bandwidth_estimators, kappa) %>%
      summarise(mise = mean(ise), med_ise=median(ise), sd_ise=sd(ise), reps= n()) %>%
      ungroup() ->
      data_mise

  }else if("lambda" %in% names(data_ise)){
    data_ise %>%
      group_by(den, n, kernel, bandwidth_estimators, lambda) %>%
      summarise(mise = mean(ise), med_ise=median(ise), sd_ise=sd(ise), reps= n()) %>%
      ungroup() ->
      data_mise
  } else {
    data_ise %>%
      group_by(den, n, kernel, bandwidth_estimators) %>%
      summarise(mise = mean(ise), med_ise=median(ise), sd_ise=sd(ise), reps= n()) %>%
      ungroup() ->
      data_mise
  }
  data_mise
}

############################################################################################
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

mise_1 <- array(mise_vec,
              dimnames = list("n_samples" =n_samples, "kernels"=c("rectangular","gaussian", "epanechnikov"),
                              "density"=c("custom density", "uniform distribution", "normal distribution")),
              dim=c(length(n_samples),length(kernels),length(dens_list)))
print(mise_1)

# -> we will set the kernel to epanechnikov and work on the largest set of samples
# the important parameter for a good estimation is the bandwidth
# above, we did set the bandwidth to h=0.1
# as you can see, the bandwidth is a much more important parameter than the kernel:
par(mfrow = c(1,1))
bandwidth_set <- list(list(0.3, "dark red"), list(0.01, "dark green"), list(0.001, "orange"))
kernel <- epanechnikov
n_samples <- 10
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
      cat("mise_bandwidths", rep)
      samples <- d[[2]](n_samples)
      kde <- kernel_density_estimator(kernel, samples, bandwidth=h[[1]])
      ise_vec <- c(ise_vec, (sum((kde$fun(grid) - dens$fun(grid))^2) * (x_lim[2] - x_lim[1])) / length(grid))
    }
    mise_vec <- c(mise_vec, mean(ise_vec))
  }
}

mise_2 <- matrix(mise_vec, nrow=c(length(bandwidth_set)),
               dimnames = list("bandwidhts"= c(0.3, 0.01, 0.001),
                               "density"=c("custom density", "uniform distribution", "normal distribution")))
print(mise_2)

#hohe reps/ns erst ab hier!

# TODO: 1. bandweitenschätzer vorstellen (mathematisch?),
# sweeter plot f?r defaultwerte auf unserer custom_density:
# TODO: main und legend mit reinwerfen (als listen?)

obj_simple_comp <- plot_comparison_objects(show_diff=FALSE, reps=200)
plot_comparison(show_diff=FALSE, reps=200, objects=obj_simple_comp)
print("plot_comparison done(616/617")
# 2. lambda/kappa wählen
# hence we estimate a upper bound for the variance, we have tuning parameters for the pco_method and goldenshluger_lepski

# first, lets take a look at goldenshluger_lepski and the kappa parameter
# in literature, they set kappa=1.2

# kappa_set <- c(1, ,1.18, 1.2, 1.22, 1.4, 1.6, 1.8, 2)
ns <- 1000
reps <- 200
kappa_set <- c(1, 1.1, 1.2, 1.3, 1.6, 2)
dens_list <- list(custom_dens=dens, dunif=Density(dunif,c(0,1)), dnorm=Density(dnorm,c(-15,15)))
sampler_list <- list(custom_sampler, runif, rnorm)
kernel_list <- list(epanechnikov=epanechnikov)
ise_kappa <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, kernels=kernel_list, kappa_set=kappa_set, ns = ns, reps = reps)
mise_kappa <- calculate_mise(ise_kappa)

print("kappa done")
# bestes kappa wird hier ... sein

# now we try to tune our lambda
# lambda_set <- c(1, 1.1, 1.2, 1.4, 1.6, 1.8, 2)
ns <- 1000
reps <- 200
lambda_set <- c(1, 1.1, 1.2, 1.4, 1.7, 2)
dens_list <- list(custom_dens=dens, dunif=Density(dunif,c(0,1)), dnorm=Density(dnorm,c(-15,15)))
sampler_list <- list(custom_sampler, runif, rnorm)
kernel_list <- list(epanechnikov=epanechnikov)
ise_lambda <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, kernels=kernel_list, lambda_set=lambda_set, ns = ns, reps = reps)
mise_lambda <- calculate_mise(ise_lambda)

print("lambda done")
# bestes lambda wird hier ... sein

# 3. Schätzer mit optimalen lambda/kappa
# nun vergleichen wir die bandweitensch?tzer untereinander
# lambda, kappa fest wie oben gew?hlt

kappa_set <- list(1.2)
lambda_set <- list(1)
reps <- 2
ns <- c(10)

# mehr densities?
# als erstes schauen wir, welcher Bandweitensch?tzer bei 1000 samples auf unseren 3 densities am besten performen w?rden
plot_object_vec_1 <- list()
for(i in seq_along(dens_list)){
  d <- dens_list[[i]]
  xlim_1 <-  c(d$support[1] - 1, d$support[2] + 1)
  plot_object_vec_1 <- c(plot_object_vec_1, list(plot_comparison_objects(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim_1[1], xlim_upper=xlim_1[2], reps=reps, ns=ns)))
}
for(i in seq_along(plot_object_vec_1)){
  obj <- plot_object_vec_1[[i]]
  par(mfrow=c(1,1))
  d <- dens_list[[i]]
  xlim_1 <-  c(d$support[1] - 1, d$support[2] + 1)
  plot_comparison(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim_1[1], xlim_upper=xlim_1[2], reps=reps, ns=ns, objects=obj)
}
ise <- compare_ise(dens_list=dens_list, dens_sampler_list=sampler_list, reps=reps,ns=ns)
mise_high_ns_comp <- calculate_mise(ise)
mise_high_ns_comp %>%
  group_by(n, bandwidth_estimators) %>%
  summarize(mean_mise=mean(mise), mean_sd_ise=mean(sd_ise))


# TODO: performance
# we will make a small performance comparison

# da nicht immer sehr viele daten vorhanden sind, betrachten wir nun, welcher Bandweitensch?tzer sich bei verschiedenen sample sizes am besten verh?lt
ns <- c(10, 50, 100, 1000)
plot_object_vec_2 <- list()
for(i in seq_along(dens_list)){
  d <- dens_list[[i]]
  xlim_2 <-  c(d$support[1] - 1, d$support[2] + 1)
  plot_object_vec_2 <- c(plot_object_vec_2, list(plot_comparison_objects(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim_2[1], xlim_upper=xlim_2[2], reps=reps, ns=ns))
)
}

# TODO: plots einzeln schöner machen
ylims <- list(c(-0.1, 2.2), c(-0.1, 1.5), c(-0.1, 0.5))
for(i in seq_along(plot_object_vec_2)){
  obj_lists <- plot_object_vec_2[[i]]
  ylims_i <- ylims[[i]]
  par(mfrow=c(1,4))
  d <- dens_list[[i]]
  if(i == 3){
    xlim_2 <-  c(-4, 4)

  }else{
    xlim_2 <-  c(d$support[1] - 1, d$support[2] + 1)

  }
  for(obj in obj_lists){
  plot_comparison(show_diff=FALSE, dens=dens_list[[i]], dens_sampler=sampler_list[[i]], xlim_lower=xlim_2[1], xlim_upper=xlim_2[2], ylim_lower=ylims_i[1], ylim_upper=ylims_i[2], reps=reps, ns=ns, objects=obj)
  }
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

save(list = ls(all.names = TRUE), file="sim_objects.rda")
#load(file="sim_objects.rda")
