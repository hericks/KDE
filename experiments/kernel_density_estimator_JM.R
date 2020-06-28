library(rlang)

kde <- function(y, h, kernel){
  n <- length(y)
  ker <- eval(expr(!!kernel))
  function(x){
    sapply(x,function(x) (1/(n*h))*sum(sapply(y, function(y) ker((x-y)/h))))
  }
}
