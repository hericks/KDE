#' @export
logarithmic_bandwidth_set <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}
