newIntegrableFunction <- function(fun, support, ..., subclass=NULL){
  stopifnot("class of fun has to be function" = class(fun) != "function")
  if(missing(support)){
    support <- find_support(fun)
  }
  stopifnot("not a valid support" = is.numeric(support) && identical(length(support), 2L))

  structure(
    list(fun, support),
    ...,
    class = c(subclass, "Integrable_function")
    )
}

IntegrableFunction <- function(fun, support){
  func <- newIntegrableFunction(fun, support)
  validatyIntegrableFunction(func)
  func
}

newDensity <- function(fun, support, ..., subclass=NULL){

  newintegrableFunction(fun,support, subclass=c(subclass,"Density"))
}

Density <- function(fun, support){
  den <- new_density(fun, support)
  validateDensity(den)
  den
}

newKernel <- function(fun, support, ..., subclass=NULL){

  newIntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}

Kernel <- function(fun, support){
  kern <- new_density(fun, support)
  validateKernel(kernel)
  kern
}


find_support <- function(fun){
  #TODO: aufräumen
  return(c(-Inf, Inf))

  x <- c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
  y <- -x
  x <- x
  y <- y

  check_x <- fun(x) != 0
  check_y <- fun(y) != 0
  if(!(any(check_x) && any(check_y))){
    return(c(-1e-6,1e-6))
  }
  non_zero <- min(c(which(check_x==TRUE), which(check_y==TRUE)))
  if(non_zero != 1){
    return(c(y[non_zero-1], x[non_zero-1]))
  }
  else{
    return(c(-Inf,Inf))
  }
}
