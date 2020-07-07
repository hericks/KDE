new_Density <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Density"))
}

Density <- function(fun, support){
  den <- new_Density(fun, support)
  validate_Density(den)
  den
}

validate_Density <- function(obj){
  if(!validate_IntegrableFunction(obj)) return(FALSE)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  neg_obj <- Vectorize(function(x) max(0, -object(x)))
  if(!isTRUE(all.equal(integrate(neg_obj, lower = lower, upper = upper)[[1]], 0))) return(FALSE)

  isTRUE(abs(integrate(object, lower = lower, upper = upper)[[1]] - 1)
         < integrate(object, lower = lower, upper = upper)[[2]])
}
