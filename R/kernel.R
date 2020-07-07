new_kernel <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}

Kernel <- function(fun, support){
  kern <- new_Kernel(fun, support)
  validate_Kernel(kern)
  kern
}

validate_Kernel <- function(obj){
  if(!validate_IntegrableFunction(obj)) return(FALSE)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]
  isTRUE(abs(integrate(object, lower = lower, upper = upper)[[1]] - 1)
         < integrate(object, lower = lower, upper = upper)[[2]])
}
