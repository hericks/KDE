evaluate <- function(fun, x) {
  UseMethod("evaluate")
}

evaluate.Integrable_function <- function(fun, x) {
  vals <- fun(x)
  vals[!is.numeric(vals) | is.na(vals)] <- 0
  vals
}

evaluate.Density <- function(fun, x) {
  vals <- fun(x)
  vals[!is.numeric(vals) | vals < 0 | is.na(vals)] <- 0
  vals
}
