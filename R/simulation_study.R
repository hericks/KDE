compare <- function(eval_points,
                    funs = list(runif),
                    bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                    ns = 50,
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

plot_comparison <- function(dens = Density(dunif, c(0, 1)),
                            dens_sampler = runif,
                            xlim_lower = -1,
                            xlim_upper = 2,
                            main = NA,
                            legend = NULL,
                            show_diff = TRUE,
                            split = TRUE,
                            bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                            reps = 4,
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
      plot(
        x_grid,
        dens_eval,
        type = "l",
        lwd = 2,
        col = 1,
        ylim = range(dens_eval) + del * c(-1, 1)
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
    plot(
      x_grid,
      dens_eval,
      type = "l",
      lwd = 2,
      col = 1,
      ylim = range(dens_eval) + del * c(-1, 1)
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

library(tidyverse)

compare_ise <- function(dens_list = list(dunif=Density(dunif, c(0, 1))),
                        dens_sampler_list = list(runif=runif),
                        bandwidth_estimators = list(cross_validation=cross_validation, goldenshluger_lepski=goldenshluger_lepski, pco_method=pco_method),
                        ns = 50,
                        kernels = list(gaussian=gaussian),
                        lambda_set = list(1),
                        kappa_set = list(1.2),
                        reps = 3,
                        num_eval_points= 300,
                        n_bandwidths = 5){

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

  print(m)
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

