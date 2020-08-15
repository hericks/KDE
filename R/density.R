#' Density
#'
#' @description
#' The S3 Class \code{Density} tries to ensure some of the properties of densities
#' and is a subclass of the \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}}
#' (see: 'Details' for exact requirements).
#'
#' @details{
#' A density function is a real valued, non-negative, integrable function, such that its integral over the real numbers equals one.
#' Density functions as \code{R} functions are required to:
#'
#'   1. be vectorised in its argument, taking a single numeric argument,
#'   returning a numerical vector of the same length only
#'
#'   2. return zero for inputs outside their support
#'
#'   3. return non-negative values for inputs inside their support
#'
#'   4. be integrated over their support using \code{integrate} without
#'   throwing an error
#'
#'   5. computate the integral properly over its support, using \code{integrate}. The result should be equal to one.
#'
#'   The functions in this package don't just take \code{R} functions satisfying
#'   these conditions, but objects of S3 class \code{IntegrableFunction} (or one
#'   of its subclasses \code{Kernel}, \code{Density}).
#'
#'   The S3 class \code{Density} exists to ensure some of the most
#'   basic properties of density functions. The class is build on lists
#'   containing two named entries \code{fun} and \code{support}.
#'
#'   * \strong{`fun`} is an \code{R} function (the represented function) taking
#'   a single numeric argument and returning a numeric vector of the same
#'   length. This function should return zero outside, and non-negative values inside of the interval given in
#'   the \code{support} entry. Using \code{integrate} over its support should evaluate to one.
#'
#'   * \strong{`support`} is a numeric vector of length 2 containing a lower-
#'   and upperbound for the support of the function stored in \code{fun} in its
#'   first and second entry respectively. In particular the values \code{-Inf}
#'   and \code{Inf} are allowed.
#'
#'   * \strong{`subdivisions`} is a integer value used for the subdivisions parameter for \code{\link[stats:integrate]{integrate}}.
#'   The function \code{fun} is needed to be integrated using \code{\link[stats:integrate]{integrate}}.
#'   Because of that, the subdivisions parameter is required to be large enough, such that \code{\link[stats:integrate]{integrate}} can work properly.
#'   The default value is set to 1000L. Be aware that too large numbers can cause long runtimes!
#'
#'   The constructor \code{Density} tries to construct a valid
#'   \code{Density} object based on the passed arguments. Returned
#'   objects are guaranteed to pass the validator [validate_Density]
#'   (\bold{Attention:} This does not guarantee the conditions in the first
#'   'Details' paragraph: see [validate_Density].).}
#'
#' @param fun an \code{R} function taking a single numeric argument and
#' returning a numeric vector of the same length: see 'Details'.
#' @param support a numerical vector of length 2 containing the lower- and
#'   upperbound in the first and second entry respectively.
#'   \code{Density} will try to find bounds on the support itself if
#'   \code{NULL} is passed.
#' @param subdivisions a integer vector of length 1, used for the subdivisions parameter for the function \code{\link[stats::integrate]{integrate}}.
#'
#' @examples
#' dens_norm <- Density(dnorm, c(-Inf, Inf))
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
#' @seealso
#' \code{\link[KDE:validate_Density]{validate_Density}}
#' \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}
#'
#' @include integrable_function.R
#'
#' @export
Density <- function(fun, support = NULL, subdivisions = 1000L){
  den <- new_Density(fun, support, subdivisions)
  validate_Density(den)
  den
}

#'Validator for S3 class \code{Density}
#'
#' @description This function serves as a validator for the S3 class
#'   \code{Density}. See 'Details' for further information and
#'   potential flaws.
#'
#' @param x an \code{R} object to validate as object of S3 class
#'   \code{[Density]}.
#'
#' @details The validator \code{validate_Density} can be used to
#'   verify objects as formally correct S3 objects of class
#'   [Density]. In particular the formal structure is ensured and it makes use of the [validate_IntegrableFunction].
#'   Additionally this function \emph{tries to} (see 'Special Attention')
#'   validate the additional conditions of valid integrable functions (as
#'   specified in the first 'Details'-paragraph of [Density]).
#'
#' @section Special Attention:
#'
#'   Like all numerical routines, \code{validate_Density} can
#'   evaluate the represented function on a finite set of points only. If the
#'   represented function returns valid results over nearly all its range, it is
#'   possible that this function misses unexpected/wrong return values. Thus,
#'   using [Density] or [validate_Density] to construct
#'   and validate objects representing integrable functions is \emph{not}
#'   sufficient to ensure the properties \[1-5\] listed in the first
#'   'Details'-paragraph of [Density], but serves more as a
#'   sanity-check.
#'
#' @seealso
#' * [Density] for more information about Density functions and the S3 class \code{Density}.
#' * [IntegrableFunction] for more information about integrable functions and the S3 class \code{IntegrableFunctions}.
#' * [validate_IntegrableFunction] for more information about the IntegrableFunction validator.
#'
#' @include integrable_function.R
#'
#'@export
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

  stopifnot("The integral of a density over its support has to be one"=
              isTRUE(all.equal(integrate(object, lower = lower, upper = upper, subdivisions = subdivisions)[[1]], 1)))
  invisible(x)
}

new_Density <- function(fun, support, subdivisions = 1000L, ..., subclass=NULL){
  new_IntegrableFunction(fun, support, subdivisions, subclass=c(subclass,"Density"))
}
