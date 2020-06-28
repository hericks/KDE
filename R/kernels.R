

rectangular <- function(x, y, h){
  u <- (x-y)/h
  return(1/2*(abs(u) <= 1))
}

triangular <- function(x, y, h){
  u <- (x-y)/h
  return((1 - abs(u)) * (abs(u) <= 1))
}

epanechnikov <- function(x, y, h){
  u <- (x-y)/h
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

biweight <- function(x, y, h){
  u <- (x-y)/h
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

triweight <- function(x, y, h){
  u <- (x-y)/h
  return(35/32 (1 - u^2)^3 * (abs(u) <= 1))
}

tricube <- function(x, y, h){
  u <- (x-y)/h
  return(70/81 (1 - abs(u^3))^3 * (abs(u) <= 1))
}

gaussian <- function(x, y, h){
  u <- (x-y)/h
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

cosine <- function(x, y, h){
  u <- (x-y)/h
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

logistic <- function(x,y,h){
  u <- (x-y)/h
  return(1/(exp(u) + 2 + exp(-u)) )
}

sigmoid_fu <- function(x,y,h){
  u <- (x-y)/h
  return(1/pi * (1/(exp(u) + exp(-u)) ) )
}


silverman <- function(x,y,h){
  u <- (x-y)/h
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}



