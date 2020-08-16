#' @include builtin_kernels.R
#' @export
compare <- function(eval_points, funs=list(runif),
                    bandwidth_estimators=list(cross_validation, goldenshluger_lepski, pco_method),
                    ns=50, kernels=list(gaussian), lambda_set=1, kappa_set=1.2, reps=400, length.out=5){
  # TODO: Argchecks necessary?
  # TODO: length.out nach Tests auf 30 setzen, reps hochsetzen
  if(!is.list(kernels)) kernels <- as.list(kernels)
  #if(!is.list(lambda_set)) lambda_set <- as.list(lambda_set)
  #if(!is.list(kappa_set)) kappa_set <- as.list(kappa_set)
  if(!is.list(ns)) ns <- as.list(ns)

  m <- length(funs) * length(ns) * length(kernels)
  m <- m*(any(sapply(bandwidth_estimators, identical, cross_validation))
       + length(kappa_set)*any(sapply(bandwidth_estimators, identical, goldenshluger_lepski))
       + length(lambda_set)*any(sapply(bandwidth_estimators, identical, pco_method)))

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
                if(any(kde$fun(eval_points) < 0)){ print(kde)}
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
                if(any(kde$fun(eval_points) < 0)){ print(kde)}

                kde$fun(eval_points)
              } )
              cnt <- cnt + 1
            }
          }else if(identical(est, cross_validation)){
            res[,,cnt] <- replicate(reps, {
              samples <- f(n)

              bandwidth <- cross_validation(k, samples, bandwidths, subdivisions)

              kde <- kernel_density_estimator(k, samples, bandwidth, subdivisions)
              if(any(kde$fun(eval_points) < 0)){ print(kde)}

              kde$fun(eval_points)
            } )
            cnt <- cnt + 1
          }
        }
      }

    }
  }
  res
}
