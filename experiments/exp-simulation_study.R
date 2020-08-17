plot_with_confidence_band <- function(x, Y, col){
  rgb <- col2rgb(col)/255
  col_alpha <- rgb(rgb[1], rgb[2], rgb[3], 0.2)
  v <- apply(Y, 1, function(x) c(mean(x), sd(x)))
  lines(x, v[1,], lwd=2, col=col)
  polygon(
    c(x, rev(x)),
    c(v[1,] + v[2,], rev(v[1, ]-v[2, ])),
    col= col_alpha,
    border=NA
  )
  lines(x, v[1, ]+v[2, ], lwd=1, col=col)
  lines(x, v[1, ]-v[2, ], lwd=1, col=col)
}

plot_comparison <- function(dens=Density(dunif,c(0,1)),
                            dens_sampler=runif,
                            xlim_lower=-1,
                            xlim_upper=2,
                            main=NA,
                            legend=NULL,
                            show_diff=TRUE,
                            split=TRUE,
                            bandwidth_estimators=list(cross_validation, goldenshluger_lepski, pco_method),
                            reps=4,
                            kappa=list(1.2),
                            lambda=list(1),
                            ...){
  dens_fun <- dens$fun
  x_grid <- seq(xlim_lower, xlim_upper, length.out=300)
  dens_eval <- dens_fun(x_grid)
  if(length(kappa) > 1 && length(lambda)>1){ stop("how is this plot supposed to look like?")}
  else if(length(kappa) > 1){
    res <- compare(x_grid, list(dens_sampler), reps=reps, bandwidth_estimators=list(goldenshluger_lepski), kappa_set=kappa, ...)
  }
  else if(length(lambda) > 1){
    res <- compare(x_grid, list(dens_sampler), reps=reps, bandwidth_estimators=list(pco_method), lambda_set=lambda, ...)
  } else {
    res <- compare(x_grid, list(dens_sampler), reps=reps, bandwidth_estimators=bandwidth_estimators, ...)
  }
  m <- dim(res)[3]
  bw_len <- length(bandwidth_estimators)
  if(isTRUE(split) && bw_len > 1){
    m_bw_len <- m / bw_len
    for(i in 1:m_bw_len){
      res_sub <- res[,,(1 + (i-1) * bw_len):(i*bw_len)]
      del <- max(apply(res_sub, c(1,3), sd))

      # plot results
      par(mar=c(0,0,1,0), ann=FALSE, xaxt="n", yaxt="n")
      plot(x_grid, dens_eval, type="l", lwd=2, col=1, ylim=range(dens_eval)+del*c(-1,1))
      grid()
      for(i in 1:bw_len){
        plot_with_confidence_band(x_grid, res_sub[,,i], col=i+1)
      }
      title(main=main)
      if(!is.null(legend)) legend("topright", legend=legend, lwd=2, col=2:(bw_len + 1))

      # second plot (difference between true f and estimation)
      if(show_diff){
        diff <- res_sub - dens_eval
        plot(c(0,1), c(0,0), type="l", lwd=2, col=1, ylim=c(-del,del))
        grid()
        for(i in 1:bw_len) plot_with_confidence_band(x_grid, diff[,,i], col=i+1)
      }
    }
  }else{
    m <- dim(res)[3]
    del <- max(apply(res,c(1,3), sd))

    # plot results
    par(mar=c(0,0,1,0), ann=FALSE, xaxt="n", yaxt="n")
    plot(x_grid, dens_eval, type="l", lwd=2, col=1, ylim=range(dens_eval)+del*c(-1,1))
    grid()
    for(i in 1:m){
      plot_with_confidence_band(x_grid, res[,,i], col=i+1)
    }
    title(main=main)
    if(!is.null(legend)) legend("topright", legend=legend, lwd=2, col=2:(m+1))

    # second plot (difference between true f and estimation)
    if(show_diff){
      diff <- res - dens_eval
      plot(c(0,1), c(0,0), type="l", lwd=2, col=1, ylim=c(-del,del))
      grid()
      for(i in 1:m) plot_with_confidence_band(x_grid, diff[,,i], col=i+1)
    }
  }
}
