# tidyverse style guide: all comments are awnsering "why" not "what" or "how"

rejection_sampling <- function(f_den, g_den, g, M) {
  force(f_den)
  force(g_den)
  force(g)
  force(M)

  # M is constant in the real numbers
  stopifnot(is.numeric(M))

  ## f_den, g_den should be probability density functions
  # f_den, g_den are nonnegative
  stopifnot(class(f_den) == "function")
  neg_f <- Vectorize(function(x) max(0, -f_den(x)))
  neg_integral_f <- integrate(neg_f, -Inf, Inf)[[1]]
  stopifnot(neg_integral_f == 0)

  stopifnot(class(g_den) == "function")
  neg_g <- Vectorize(function(x) max(0, -g_den(x)))
  neg_integral_g <- integrate(neg_g, -Inf, Inf)[[1]]
  stopifnot(neg_integral_g == 0)

  #f_den, g_den are normalized
  pos_f <- Vectorize(function(x) max(0, f_den(x)))
  pos_integral_f <- integrate(pos_f, -Inf, Inf)[[1]]
  stopifnot(pos_integral_f == 1)

  pos_g <- Vectorize(function(x) max(0, g_den(x)))
  pos_integral_g <- integrate(pos_g, -Inf, Inf)[[1]]
  stopifnot(pos_integral_g == 1)

  # relation between f_den,M and g_den that has to be satisfied
  temp <- runif(1e6, -1e12, 1e12)
  for (v in temp) {
    if (f_den(v) > M*g_den(v)) {
      return(FALSE)
    }
  }

  function(n) {
    # n should be a nonnegative integer
    stopifnot(n%%1 == 0)
    stopifnot(n >= 0)

    u <- numeric(M*n)
    samples <- numeric(M*n)
    accepted <- NULL

    while(length(accepted) < n) {
      u <- runif(M*n)
      samples <- g(M*n)

      accepted <- c(accepted,samples[u*M*g_den(samples) < f_den(samples)])
    }

    accepted[1:n]
  }
}
