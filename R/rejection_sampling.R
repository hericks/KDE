rejection_sampling <- function(f_den, g_den, g, M) {
  # TODO: Check sufficient conditions on arguments

  force(f_den)
  force(g_den)
  force(g)
  force(M)

  function(n) {
    # TODO: Check sufficient conditions on arguments

    u <- numeric(M*n)
    samples <- numeric(M*n)
    accepted <- NULL

    while(length(accepted) < n) {
      u <- runif(2*n)
      samples <- g(2*n)

      accepted <- samples[u*M*g_den(samples) < f_den(samples)]
    }

    accepted[1:n]
  }
}
