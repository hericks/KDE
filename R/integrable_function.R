#' Integrable Functions
#'
#' @description Many functions of the \code{KDE} package work with densities and
#'   kernels, which are integrable in the mathematical sense. See 'Details' for
#'   exact requirements. The S3 class \code{IntegrableFunction} tries to ensure
#'   some of the properties of integrable functions and serves as superclass for
#'   the more specific S3 classes \code{Density} and \code{Kernel}.
#'
#' @param fun a \code{R} function taking a single numeric argument and returning
#'   a numeric vector of the same length. See 'Details' for further
#'   requirements.
#' @param support numerical vector of length 2; the lower- and upper bound of
#'   the compact support in the first and second entry respectively. In
#'   particular non-finite values are prohibited. \code{IntegrableFunction} will
#'   try to find bounds on the support itself if \code{NULL} is passed.
#' @param subdivisions positive numeric scalar; the subdivisions parameter for
#'   the function \code{\link{integrate_primitive}}.
#' @param ... additional parameters to keep fixed during the evaluation of
#'   \code{fun}.
#'
#' @details Integrable functions as \code{R} functions are required to
#'
#'   1. be vectorised in its argument, taking a single numeric argument,
#'   returning a numerical vector of the same length only,
#'
#'   2. return zero for inputs outside their compact support,
#'
#'   3. can be integrated over their support using \code{integrate_primitive}
#'   and the given number of subdivisions (the relative error converges).
#'
#'   Notice that a compact support may sound like a strong restriction, but
#'   since every integrable function is near zero outside of a compact set this
#'   is computatianlly always given for integrable functions.
#'
#'   The functions in this package don't just take \code{R} functions satisfying
#'   these conditions, but objects of S3 class \code{IntegrableFunction} (or one
#'   of its subclasses \code{Kernel}, \code{Density}).
#'
#'   The S3 class \code{IntegrableFunction} exists to ensure some of the most
#'   basic properties of integrable functions. The class is build on lists
#'   containing three named entries \code{fun}, \code{support} and \code{subdivisions}
#'
#'   * \strong{`fun`} is a \code{R} function (the represented function) taking a
#'   single numeric argument (additional to the fixed arguments in \code{...})
#'   and returning a numeric vector of the same length. This function should
#'   return near zero outside of the interval given in the \code{support} entry.
#'
#'   * \strong{`support`} is a numeric vector of length 2 containing a lower-
#'   and upperbound for the support of the function stored in \code{fun} in its
#'   first and second entry respectively. In particular the values \code{-Inf}
#'   and \code{Inf} are allowed.
#'
#'   * \strong{`subdivisions`} is a positive numeric scalar used as the
#'   subdivisions parameter for \code{\link{integrate_primitive}}. The function
#'   \code{fun} should be integrable (the relative error converges). Therefore
#'   the subdivisions parameter should to be large enough, such that
#'   \code{integrate_primitive} yields a sufficiently accurate result.
#'
#'   The constructor \code{IntegrableFunction} tries to construct a valid
#'   \code{IntegrableFunction} object based on the passed arguments. Returned
#'   objects are guaranteed to pass the validator \code{validate_IntegrableFunction}.
#'
#'   \bold{Attention:} This does not guarantee the conditions in the first
#'   'Details' paragraph. See \code{\link{validate_IntegrableFunction}} for further
#'   information.
#'
#' @seealso \code{\link{integrate_primitive}} for the integration method used,
#'   \code{\link{Kernel}}/\code{\link{Density}} for more information about
#'   kernels/densities.
#'
#' @include integrate_primitive.R
#'
#' @export
IntegrableFunction <- function(fun, support=NULL, subdivisions=1000L, ...){
  extra_args <- list(...)
  new_fun <- function(x) do.call(fun, c(list(x), extra_args))
  obj <- new_IntegrableFunction(new_fun, support, subdivisions)
  validate_IntegrableFunction(obj)
  obj
}

#' Validator for S3 class \code{IntegrableFunction}
#'
#' @description This function serves as a validator for the S3 class
#'   \code{IntegrableFunction}. See 'Details' for further information and
#'   potential flaws.
#'
#' @param x an \code{R} object to validate as object of S3 class
#'   \code{IntegrableFunction}.
#'
#' @details The validator \code{validate_IntegrableFunction} can be used to
#'   verify objects as formally correct S3 objects of class
#'   \code{\link{IntegrableFunction}}. In particular the formal structure is ensured.
#'   Additionally this function \emph{tries to} (see 'Special Attention')
#'   validate the additional conditions of valid integrable functions as
#'   specified in the first 'Details'-paragraph of \code{IntegrableFunction}.
#'
#' @section Special Attention:
#'
#'   Like all numerical routines, \code{validate_IntegrableFunction} can
#'   evaluate the represented function on a finite set of points. If the
#'   represented function returns valid results over nearly all its range, it is
#'   possible that the validator misses unexpected/wrong return values. Thus,
#'   using \code{IntegrableFunction} or \code{validate_IntegrableFunction} to construct
#'   and validate objects representing integrable functions is \emph{not}
#'   sufficient to ensure the properties listed in the first
#'   'Details'-paragraph of \code{IntegrableFunction}, but serves more as a
#'   sanity-check.
#'
#' @seealso
#' \code{\link{IntegrableFunction}} for more information about integrable functions and the S3 class \code{IntegrableFunction}.
#'
#' @include integrate_primitive.R
#'
#' @export
validate_IntegrableFunction <- function(x){
  # to test the basic structure
  stopifnot("Object must inherit 'IntegrableFunction'"=inherits(x, "IntegrableFunction"))
  stopifnot("Object must be a list"=is.list(x))
  stopifnot("Object must contain entry 'fun'"="fun" %in% names(x))
  stopifnot("Object must contain entry 'support'"="support" %in% names(x))
  stopifnot("Object must contain entry 'subdivisions'"="subdivisions" %in% names(x))

  # to test the structure of 'fun' entry
  fun <- x$fun
  stopifnot("Entry 'fun' must be a function"=is.function(fun))

  # to test the structure of 'support' entry
  support <- x$support
  stopifnot("Entry 'support' must be numeric"=is.numeric(support))
  stopifnot("Entry 'support' must be of length 2"=identical(length(support), 2L))
  stopifnot("Entry 'support' must contain lower- before upperbound"=support[1] < support[2])
  stopifnot("Entry 'support' must only contain finite values"=is.finite(support[1]) & is.finite(support[2]))

  offsets <- 10**(-1:10)
  testing_values <- c(fun(max(support[1],-1e10) - offsets), fun(min(support[2], 1e10) + offsets))
  stopifnot("Entry 'fun' has to be zero outside of support"=all(sapply(testing_values, function(x) isTRUE(all.equal(x, 0)))))

  # to test the structure of the 'subdivisions' entry
  subdivisions <- x$subdivisions
  stopifnot("Entry 'subdivisions' must be numeric"=is.numeric(subdivisions))
  stopifnot("Entry 'subdivisions' must be of length 1"=length(subdivisions) == 1)
  stopifnot("Entry 'subdivisions' must be positive"=subdivisions > 0)

  # to test the the integrability over the given support using the given number of subdivisions
  integration_value <- integrate_primitive(function(x) abs(fun(x)), support[1], support[2], subdivisions = subdivisions)$value

  stopifnot("The integral of the absolute function must integrate to a finite value"=is.finite(integration_value))
  invisible(x)
}

# new_IntegrableFunction is only called by the public IntegrableFunction
# constructor and private subclass constructors. Always followed by a call of
# validate_IntegrableFunction (maybe inside a subclass validator). Further
# checking is done by validate_IntegrableFunction.
new_IntegrableFunction <- function(fun, support=NULL, subdivisions=1000L, ..., subclass=NULL){
  stopifnot("class of fun has to be function"=is.function(fun))

  if(is.null(support)){
    support <- find_support(fun)
  }

  structure(
    list("fun"=fun, "support"=support, "subdivisions"=subdivisions),
    ...,
    class=c(subclass, "IntegrableFunction")
  )
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

#' Print objects of S3 class \code{IntegrableFunction}
#'
#' @param x object of S3 class \code{IntegrableFunction}; the object to print.
#' @param class_prefix optional character vector of length 1; string to replace
#'   \code{IntegrableFunction} in formatted output.
#' @param ... unused argument; used for compatibility with generic print.
#'
#' @export
print.IntegrableFunction <- function(x, class_prefix=NULL, ...) {
  extra_args <- with(environment(x$fun), extra_args)
  prefixes <- sapply(names(extra_args), function(name) {if (name == "") "" else paste0(name, " = ")}, USE.NAMES = FALSE)
  parts <- sapply(seq_len(length(prefixes)), function(i) paste0(prefixes[i], extra_args[[i]]))
  eval_expr <- paste0("(", paste(c("x", parts), collapse=", "), ")")

  if (is.null(class_prefix)) cat("IntegrableFunction", sep="\n")
  else cat(class_prefix[1], sep="\n")

  fun_lines <- deparse(with(environment(x$fun), fun))
  collapse <- ifelse(length(fun_lines) == 2, "", "\n")
  cat(paste(fun_lines, collapse=collapse))
  cat(paste0("\nevaluated at ", eval_expr, ", support: [", x$support[1], ",", x$support[2], "], subdivisions: ", x$subdivisions))
}
