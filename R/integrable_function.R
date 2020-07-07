newIntegrableFunction <- function(fun, support, ..., subclass=NULL){
  stopifnot("class of fun has to be function" = (class(fun) == "function"))
  if(missing(support)){
    if(isFALSE(tryCatch({find_borders(object,center)}, error=function(e) FALSE))){
      stop("cannot find a support")
    }
    support <- find_support(fun)
  }
  stopifnot("not a valid support" = is.numeric(support) && identical(length(support), 2L))

  structure(
    list("fun" = fun, "support" = support),
    ...,
    class = c(subclass, "Integrable_function")
    )
}

IntegrableFunction <- function(fun, support){
  func <- newIntegrableFunction(fun, support)
  validateIntegrableFunction(func)
  func
}

newDensity <- function(fun, support, ..., subclass=NULL){

  newintegrableFunction(fun,support, subclass=c(subclass,"Density"))
}

Density <- function(fun, support){
  den <- newDensity(fun, support)
  validateDensity(den)
  den
}

newKernel <- function(fun, support, ..., subclass=NULL){

  newIntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}

Kernel <- function(fun, support){
  kern <- newKernel(fun, support)
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

validateIntegrableFunction <- function(obj){
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  # neg/pos part of the kernel function to be tested
  pos_object <- Vectorize(function(x) max(0, object(x)))
  neg_object <- Vectorize(function(x) max(0, -object(x)))

  # is_kernel returns a boolean and is not supposed to throw an error
  if(isFALSE(tryCatch({
    pos_int <- integrate(pos_object, lower=lower, upper=upper)
    neg_int <- integrate(neg_object, lower=lower, upper=upper)
  }, error=function(e) FALSE))) {
    return(FALSE)
  }
  pos_integral <- integrate(pos_object, lower, upper)[[1]]
  neg_integral <- integrate(neg_object, lower, upper)[[1]]
  if(is.infinite(pos_integral) && is.infinite(pos_integral)){ return(FALSE)}

  return(TRUE)
}

validateDensity <- function(obj){
  if(!validateIntegrableFunction(obj)) return(FALSE)
  object <- obj[[1]]
  lower <- obj[[2]][1]
  upper <- obj[[2]][2]

  neg_obj <- Vectorize(function(x) max(0, -obj(x)))
  if(!isTRUE(all.equal(integrate(neg_ob, lower = lower, upper = upper)[[1]], 0))) return(FALSE)

  isTRUE(abs(integrate(object, lower = lower, upper = upper)[[1]] - 1)
         < integrate(object, lower = lower, upper = upper)[[2]])
}

validateKernel <- function(obj){
  if(!validateIntegrableFunction(obj)) return(FALSE)
  object <- obj[[1]]
  lower <- obj[[2]][1]
  upper <- obj[[2]][2]
  isTRUE(abs(integrate(object, lower = lower, upper = upper)[[1]] - 1)
         < integrate(object, lower = lower, upper = upper)[[2]])
}

ker <- Kernel(rectangular, c(-Inf,Inf))
