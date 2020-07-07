new_Density <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Density"))
}

Density <- function(fun, support){
  den <- new_Density(fun, support)
  validate_Density(den)
  den
}

validate_Density <- function(obj){
  if(!inherits(obj, "Density")) return(FALSE)
  if(!validate_IntegrableFunction(obj)) return(FALSE)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  isTRUE(abs(integrate(evaluate_safe(obj), lower = lower, upper = upper)[[1]] - 1)
         < integrate(evaluate_safe(obj), lower = lower, upper = upper)[[2]])
}
