#' Densities
#'
#' @description The S3 class \code{Density} tries to ensure some of the
#' properties of densities and is a subclass \code{\link{IntegrableFunction}}
#' (see: 'Details' for exact requirements).
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
#' @details A density function is a real valued, non-negative, integrable
#'   function, such that its integral over the real numbers equals one. Density
#'   functions as \code{R} functions are required to
#'
#'   1. be vectorised in its argument, taking a single numeric argument,
#'   returning a numerical vector of the same length only,
#'
#'   2. return zero for inputs outside their compact support,
#'
#'   3. return non-negative values for inputs inside their support,
#'
#'   4. can be integrated over their support using \code{integrate_primitive}
#'   and the given number of subdivisions (the relative error converges),
#'
#'   5. yield an integral of nearly 1 (absolute error < 1%).
#'
#'   See the 'Details' section of \code{\link{IntegrableFunction}} for comments
#'   on the restrictiveness of compact supports.
#'
#'   The S3 class \code{Density} exists to ensure some of the most basic
#'   properties of density functions. The class is build on
#'   \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}} and inherits its
#'   structure.
#'
#'   The constructor \code{Density} tries to construct a valid \code{Density}
#'   object based on the passed arguments. Returned objects are guaranteed to
#'   pass the validator \code{\link{validate_Density}}.
#'
#'   \bold{Attention:} This does not guarantee the conditions in the first
#'   'Details' paragraph: see \code{validate_Density}.
#'
#' @examples
#' dens_norm <- Density(dnorm, c(-15, 15))
#' dens_unif <- Density(dunif)
#' x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, dens_norm$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Densities", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x,dens_unif$fun(x), col="red")
#' legend("topright",
#'        legend=c("dens_norm", "dens_unif"),
#'        col=c("black","red"), lty=1, cex=0.8)
#'
#' @seealso \code{\link{validate_Density}} for the corresponding validator,
#'   \code{\link{IntegrableFunction}} for more information about the superclass.
#'
#' @include integrable_function.R
#'
#' @export
Density <- function(fun, support = NULL, subdivisions = 1000L, ...){
  extra_args <- list(...)
  new_fun <- function(x) do.call(fun, c(list(x), extra_args))
  den <- new_Density(new_fun, support, subdivisions)
  validate_Density(den)
  den
}

#' Validator for the S3 class \code{Density}
#'
#' @description This function serves as a validator for the S3 class
#'   \code{\link{Density}}. See 'Details' for further information and potential
#'   flaws.
#'
#' @param x an \code{R} object to validate as object of S3 class \code{Density}.
#'
#' @details The validator \code{validate_Density} can be used to verify objects
#'   as formally correct S3 objects of class \code{\link{Density}}. In
#'   particular the formal structure is ensured and
#'   \code{\link{validate_IntegrableFunction}} is called internally.
#'   Additionally this function \emph{tries to} (see 'Special Attention')
#'   validate the additional conditions of valid density functions (as specified
#'   in the first 'Details'-paragraph of \code{Density}).
#'
#' @section Special Attention:
#'
#'   Like all numerical routines, \code{validate_Density} can evaluate the
#'   represented function on a finite set of points only. If the represented
#'   function returns valid results over nearly all its range, it is possible
#'   that this function misses unexpected/wrong return values. Thus, using
#'   [Density] or [validate_Density] to construct and validate objects
#'   representing integrable functions is \emph{not} sufficient to ensure the
#'   properties \[1-5\] listed in the first 'Details'-paragraph of [Density],
#'   but serves more as a sanity-check.
#'
#' @seealso \code{\link{Density}} for more information about density functions and
#' the S3 class \code{Density}, \code{\link{IntegrableFunction}} for more
#' information about the superclass \code{IntegrableFunction}.
#'
#' @include integrable_function.R
#'
#' @export
validate_Density <- function(x){
  stopifnot("object has to be of class Density"=inherits(x, "Density"))
  validate_IntegrableFunction(x)
  object <- x$fun
  lower <- x$support[1]
  upper <- x$support[2]
  subdivisions <- x$subdivisions

  testing_points <- seq(max(lower, -1e10), min(upper, 1e10), length.out=1e4)
  # prevent double checking with validate_integrable function
  stopifnot("density functions are non-negative" = all(object(testing_points) >= 0))

  I <- integrate_primitive(object, lower = lower, upper = upper, subdivisions = subdivisions)$value
  if (isTRUE(abs(I-1) > 0.01))
    stop(paste0("The integral of a kernel over its support has to be one (is ", I, ")"))

  invisible(x)
}

new_Density <- function(fun, support, subdivisions = 1000L, subclass=NULL){
  new_IntegrableFunction(fun, support, subdivisions, subclass=c(subclass,"Density"))
}

#' Print objects of S3 class \code{Density}
#'
#' @param x object of S3 class \code{Density}; the object to print.
#' @param ... unused argument; used for compatibility with generic print.
#'
#' @export
print.Density <- function(x, ...) {
  print.IntegrableFunction(x, "Density")
}
