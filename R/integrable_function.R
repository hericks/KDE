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

find_support <- function(fun) {
  testing_points <- c(-10**(10:-10), 10**(-10:10))
  testing_values <- fun(testing_points)
  testing_values[!is.numeric(testing_values)] <- 0

  non_zero_indices <- which(testing_values != 0)

  if (length(non_zero_indices) == 0L) return(c(-Inf, Inf))

  lower_index <- min(non_zero_indices) - 1
  upper_index <- max(non_zero_indices) + 1

  lower_bound <- ifelse(lower_index < 1, -Inf, testing_points[lower_index])
  upper_bound <- ifelse(upper_index > length(testing_points), Inf, testing_points[upper_index])

  c(lower_bound, upper_bound)
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
