new_IntegrableFunction <- function(fun, support, ..., subclass=NULL){
  stopifnot("class of fun has to be function" = (class(fun) == "function"))

  if(missing(support)){
      support <- find_support(fun)
  }

  stopifnot("support has to be numeric"=is.numeric(support))
  stopifnot("support has to be of length 2"=identical(length(support), 2L))
  stopifnot("support has to consist of lower- and upperbound in increasing order"=(support[1] < support[2]))

  structure(
    list("fun"=fun, "support"=support),
    ...,
    class=c(subclass, "IntegrableFunction")
    )
}

#' TODO: documentation
#' @export
IntegrableFunction <- function(fun, support){
  func <- new_IntegrableFunction(fun, support)
  validate_IntegrableFunction(func)
  func
}

find_support <- function(fun) {
  testing_points <- c(-10**(10:-10), 10**(-10:10))
  testing_values <- fun(testing_points)

  stopifnot("fun has to return numerical values"=is.numeric(testing_values))
  stopifnot("fun hat to be vectorised"=identical(length(testing_values), length(testing_points)))

  # can't use isFALSE, since all.equal return value is TRUE or a vector of mode "character"
  non_zero_indices <- which(sapply(testing_values, function(x) !isTRUE(all.equal(x, 0))))

  if (length(non_zero_indices) == 0L) return(c(-Inf, Inf))

  lower_index <- min(non_zero_indices) - 1
  upper_index <- max(non_zero_indices) + 1

  lower_bound <- ifelse(lower_index < 1, -Inf, testing_points[lower_index])
  upper_bound <- ifelse(upper_index > length(testing_points), Inf, testing_points[upper_index])

  c(lower_bound, upper_bound)
}

#' TODO: documentation
#'
#'@export
validate_IntegrableFunction <- function(obj){
  if(!inherits(obj, "IntegrableFunction")) return(FALSE)

  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  # neg/pos part of the kernel function to be tested
  pos_object <- Vectorize(function(x) max(0, object(x)))
  neg_object <- Vectorize(function(x) max(0, -object(x)))

  # is_kernel returns a boolean and is not supposed to throw an error
  if(isFALSE(tryCatch({
    pos_integral <- integrate(pos_object, lower=lower, upper=upper)[[1]]
    neg_integral <- integrate(neg_object, lower=lower, upper=upper)[[1]]
  }, error=function(e) FALSE))) {
    return(FALSE)
  }

  if(is.infinite(pos_integral) && is.infinite(pos_integral)){ return(FALSE)}

  return(TRUE)
}

