library(tidyverse)

compare_ise <- function(dens_list = list(dunif=Density(dunif, c(0, 1))),
                        dens_sampler_list = list(runif),
                        bandwidth_estimators = list(cross_validation=cross_validation, goldenshluger_lepski=goldenshluger_lepski, pco_method=pco_method),
                        ns = 50,
                        kernels = list(gaussian=gaussian),
                        lambda_set = list(1),
                        kappa_set = list(1.2),
                        reps = 3,
                        num_eval_points= 30,
                        length.out_bandwidths = 5){

  eval_points_set <- list()
  for(d in dens_list){
    eval_points_range <-  c(d$support[1] - 1, d$support[2] + 1)
    eval_points <- seq(from=eval_points_range[1],to=eval_points_range[2], length.out=num_eval_points)
    eval_points <- list(eval_points)
    eval_points_set <- c(eval_points_set, eval_points)
  }
  time <-system.time(res <- compare(eval_points=eval_points_set, funs=dens_sampler_list, ns=ns, kernels=kernels,
                                   bandwidth_estimators=bandwidth_estimators,
                                   lambda_set=lambda_set, kappa_set=kappa_set, reps=reps, length.out=length.out_bandwidths))
  print(time)

  #f_true <- sapply(dens_list, function(d) d$fun(x_grid))
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
  ise <-apply(diff^2,c(2, 3), mean)

  if (length(kappa_set) > 1 &&
      length(lambda_set) > 1) {
    stop("how is this supposed to look like?")
  }
  else if (length(kappa_set) > 1) {
    opts <-as_tibble(expand.grid(kappa_set=kappa_set, bandwidth_estimators=bandwidth_estimators, kernel=names(kernels), n=ns, dens=names(dens_list)))
  }
  else if (length(lambda_set) > 1) {
    opts <-as_tibble(expand.grid(lambda_set=lambda_set, bandwidth_estimators=bandwidth_estimators, kernel=names(kernels), n=ns, dens=names(dens_list)))

  } else {
    opts <-as_tibble(expand.grid(bandwidth_estimators=names(bandwidth_estimators), kernel=names(kernels), n=ns, dens=names(dens_list)))
  }
  add_column(opts[rep(1:nrow(opts), each=reps), ], ise=as.vector(ise))
}
dens_list <- list(
  uniform_distri=Density(dunif, c(0,1)),
  normal_distri=Density(dnorm, c(-15,15))

)

sampler_list <- list(
  uniform_distri=runif,
  normal_distri=rnorm
)

data_ise <- compare_ise(dens_list, sampler_list)

#bandwidth_estimators=bandwidth_estimators, kernel=names(kernels), n=ns, dens=names(dens_list)

overall_mise <- function(data_ise){
  if("kappa_set" %in% names(data_ise)){

  }else if("lambda_set" %in% names(data_ise)){

  } else {

  }
}
