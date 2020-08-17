#' @include builtin_kernels.R
#' @export
compare <- function(eval_points,
                    funs = list(runif),
                    bandwidth_estimators = list(cross_validation, goldenshluger_lepski, pco_method),
                    ns = 50,
                    kernels = list(gaussian),
                    lambda_set = list(1),
                    kappa_set = list(1.2),
                    reps = 400,
                    length.out = 5) {
  if (!is.list(kernels))
    kernels <- as.list(kernels)
  if (!is.list(ns))
    ns <- as.list(ns)

  if(is.list(eval_points)){
    if(!(length(unique(sapply(eval_points, length))) == 1)) stop("same length for all eval_points vectors are needed")
    if(length(eval_points) != length(funs) && length(eval_points) > 1) stop("eval_point set has to be of same length as funs")
  }


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
  if (length(kappa) > 1 &&
      length(lambda) > 1) {
    stop("how is this plot supposed to look like?")
  }
  else if (length(kappa) > 1) {
    res <-
      compare(
        x_grid,
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
        x_grid,
        list(dens_sampler),
        reps = reps,
        bandwidth_estimators = list(pco_method),
        lambda_set = lambda,
        ...
      )
  } else {
    res <-
      compare(
        x_grid,
        list(dens_sampler),
        reps = reps,
        bandwidth_estimators = bandwidth_estimators,
        ...
      )
  }
  m <- dim(res)[3]
  bw_len <- length(bandwidth_estimators)
  if (isTRUE(split) && bw_len > 1) {
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
