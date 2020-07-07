evaluate <- function(fun, x) {
  UseMethod("evaluate")
}

evaluate.IntegrableFunction <- function(obj, x) {
  fun <- obj$fun
  vals <- fun(x)
  vals[!is.numeric(vals) | is.na(vals)] <- 0
  vals
}

evaluate.Density <- function(obj, x) {
  fun <- obj$fun
  vals <- fun(x)
  vals[!is.numeric(vals) | vals < 0 | is.na(vals)] <- 0
  vals
}

# get "safe" function out of object
evaluate_safe <- function(fun, x) {
  UseMethod("evaluate_safe")
}

evaluate_safe.IntegrableFunction <- function(obj){
  fun <- obj$fun
  function(x){
    vals <- fun(x)
    vals[!is.numeric(vals) | is.na(vals)] <- 0
    vals
  }
}

evaluate_safe.Density <- function(obj){
  fun <- obj$fun
  function(x){
    vals <- fun(x)
    vals[!is.numeric(vals) | vals < 0 | is.na(vals)] <- 0
    vals
  }
}
