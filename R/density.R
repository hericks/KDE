new_Density <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Density"))
}


#'
#'
#'@export
Density <- function(fun, support){
  den <- new_Density(fun, support)
  validate_Density(den)
  den
}

validate_Density <- function(obj){
  stopifnot("object has to be of class Density"=inherits(obj, "Density"))
  validate_IntegrableFunction(obj)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  stopifnot("The integral of a kernel over its support has to be one"=
              (abs(integrate(evaluate_safe(obj), lower = lower, upper = upper)[[1]] - 1)
               < integrate(evaluate_safe(obj), lower = lower, upper = upper)[[2]]))
  invisible(obj)
}
