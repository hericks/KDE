#' Construct a Logarithmic Bandwidth Set
#'
#' @description \code{logarithmic_bandwidth_set} returns a logarithmically
#'   scaled sequence of \code{length.out} elements between \code{from} and
#'   \code{to}.
#'
#' @param from strictly positive numeric scalar; the starting value of the
#'   sequence.
#' @param to strictly positive numeric scalar; the end value of the sequence.
#' @param length.out non-negative numeric scalar; desired length of the sequence
#'   (will be rounded up if fractional).
#'
#' @details Numerical inputs should all be \link[base:is.finite]{finite} (that is,
#'   not infinite, \code{NaN} or \code{NA}).
#'
#' @return A logarithmically spaced sequence of \code{length.out} elements
#'   between \code{from} and \code{to}.
#'
#' @seealso \code{\link{seq}} for linearly scaled sequences.
#'
#' @export
logarithmic_bandwidth_set <- function(from, to, length.out) {
  # Conditions on 'from' argument
  stopifnot("from must be numeric"=is.numeric(from))
  stopifnot("from must be of length 1"=isTRUE(all.equal(length(from), 1)))
  stopifnot("from must be strictly positive"=isTRUE(from > 0))
  stopifnot("from must be finite"=is.finite(from))

  # Conditions on 'to' argument
  stopifnot("to must be numeric"=is.numeric(to))
  stopifnot("to must be of length 1"=isTRUE(all.equal(length(to), 1)))
  stopifnot("to must be strictly positive"=isTRUE(to > 0))
  stopifnot("to must be finite"=is.finite(from))

  # Conditions on 'length.out' argument
  stopifnot("length.out must be numeric"=is.numeric(length.out))
  stopifnot("length.out must be of length 1"=isTRUE(all.equal(length(length.out), 1)))
  stopifnot("length.out must be strictly positive"=isTRUE(length.out > 0))
  stopifnot("length.out must be finite"=is.finite(length.out))

  exp(seq(log(from), log(to), length.out = ceiling(length.out)))
}
