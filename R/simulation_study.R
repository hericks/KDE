compare <- function(eval_points, funs=list(runif), ns=50, kernels=gaussian, lambda_set=1, kappa_set=1.2, reps=400, length.out=5){
  # TODO: Argchecks necessary?
  # TODO: length.out nach Tests auf 30 setzen
  if(!is.list(kernels)) kernels <- list(kernels)
  if(!is.list(bandwidth_estimators)) bandwidth_estimators <- list(bandwidth_estimators)
  if(!is.list(lambda_set)) lambda_set <- list(lambda_set)
  if(!is.list(kappa_set)) kappa_set <- list(kappa_set)


  bandwidth_estimators <- list(cv=cross_validation, gl=goldenshluger_lepski, pco=pco_method)

  m <- length(funs) * length(ns) * length(kernels)
  m <- m * length(lambda_set) + m * length(kappa_set) + m
  res <- array(NA, dim=c(length(eval_points), reps, m))
  cnt <- 1
  subdivisions <- 1000L

  for(f in funs){
    for(n in ns){
      bandwidths <- logarithmic_bandwidth_set(from=1/n, to=1, length.out=length.out)
      for(k in kernels){
        for(est in bandwidth_estimators){
          if(identical(est, goldenshluger_lepski)){
            for(kappa in kappa_set){
              res[,,cnt] <- replicate(reps, {
                samples <- f(n)

                bandwidth <- goldenshluger_lepski(k, samples, bandwidths, kappa, subdivisions)

                kde <- kernel_density_estimator(k, samples, bandwidth, subdivisions)
                kde$fun(eval_points)
              } )
              cnt <- cnt + 1
            }
          }else if(identical(est, pco_method)){
            for(lambda in lambda_set){
              res[,,cnt] <- replicate(reps, {
                samples <- f(n)

                bandwidth <- pco_method(k,samples, bandwidths, lambda, subdivisions)

                kde <- kernel_density_estimator(k, samples, bandwidth, subdivisions)
                kde$fun(eval_points)
              } )
              cnt <- cnt + 1
            }
          }else {
            for(lambda in lambda_set){
              res[,,cnt] <- replicate(reps, {
                samples <- f(n)

                bandwidth <- cross_validation(k,samples, bandwidths, subdivisions)

                kde <- kernel_density_estimator(k, samples, bandwidth, subdivisions)
                kde$fun(eval_points)
              } )
              cnt <- cnt + 1
            }
          }
        }
      }

    }
  }
  res
}
