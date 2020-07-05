


checkKernelConditions <- function(u){
  stopifnot("the argument of a kernel has to be numeric"=is.numeric(u))
}

rectangular <- function(u){
  checkKernelConditions(u)
  return(1/2*(abs(u) <= 1))
}

triangular <- function(u){
  checkKernelConditions(u)
  return((1 - abs(u)) * (abs(u) <= 1))
}

epanechnikov <- function(u){
  checkKernelConditions(u)
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

biweight <- function(u){
  checkKernelConditions(u)
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

triweight <- function(u){
  checkKernelConditions(u)
  return((35/32 * (1 - u^2)^3) * (abs(u) <= 1))
}

tricube <- function(u){
  checkKernelConditions(u)
  return(70/81 * (1 - abs(u^3))^3 * (abs(u) <= 1))
}

gaussian <- function(u){
  checkKernelConditions(u)
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

cosine <- function(u){
  checkKernelConditions(u)
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

logistic <- function(u){
  checkKernelConditions(u)
  return(1/(exp(u) + 2 + exp(-u)) )
}

sigmoidFunction <- function(u){
  checkKernelConditions(u)
  return(2/pi * (1/(exp(u) + exp(-u)) ) )
}

silverman <- function(u){
  checkKernelConditions(u)
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}



