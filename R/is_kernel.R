is_kernel <- function(kernel) {
  if (class(kernel) != "function") return(FALSE)

  pos_kernel <- Vectorize(function(x) max(0, kernel(x)))
  neg_kernel <- Vectorize(function(x) max(0, -kernel(x)))

  pos_integral <- integrate(pos_kernel, -Inf, Inf)[[1]]
  neg_integral <- integrate(neg_kernel, -Inf, Inf)[[1]]

  if (is.infinite(pos_integral) && is.infinite(pos_integral)) return(FALSE)

  isTRUE(all.equal(pos_integral - neg_integral, 1))
}
