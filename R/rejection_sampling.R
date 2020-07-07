#'Rejection Sampling (I.e., Accept-reject Method)
#'
#'\code{rejection_sampling} returns a function that draws samples from a
#'difficult probability density function.
#'
#'@param f_den normalized probability density function, from which the desired
#'  samples should be drawn.
#'@param g_den normalized probability density function that is used as a
#'  help-function to draw samples from, instead of drawing samples from
#'  \code{f_den}.
#'@param g random number generator function that draws samples from the densitiy
#'  function \code{g_den}.
#'@param M real number, which satisfies the condition that \code{f(x)} is less
#'  or equal than \code{M*g(x)} for all real numbers x.
#'@details the algorithm first draws samples from the more "well-behaved"
#'  density function \code{g_den} and then accepts/rejects these draws,
#'  according to whether they are likely within the proposed denisty function
#'  \code{f_den} or not. It is the most efficient if \code{M} is chosen as small
#'  as possible. That means that the shape of g_den is very close to f_den and
#'  as a result it is more likely that the samples drawn from \code{g_den} get
#'  accepted.
#' @examples f_den <- function(x) {
#' ret <- 1 + sin(2*pi*x)
#' ret[x < 0 | 1 < x] <- 0
#' ret}
#'custom_sampler <- rejection_sampling(f_den, dunif, runif, 2)
#'x <- seq(-0.5, 1.5, by=0.01)
#'y <- f_den(x)
#'plot(x, y, type="l", main="Custom density: 1 + sin(2*pi*x)", ylab="density")
#'n <- 65
#'points(custom_sampler(65), rep(1, n), col="red", cex=0.5)
#'
#'@return returning a function with argument \code{n} that draws \code{n}
#'  samples from the density function \code{f_den}.
#'@export


rejection_sampling <- function(f_den, g_den, g, M) {
  force(f_den)
  force(g_den)
  force(g)
  force(M)

  # M is constant in the real numbers
  stopifnot(is.numeric(M))

  # f_den, g_den should be probability density functions
  stopifnot(is_density(g_den) == TRUE)
  stopifnot(is_density(g_den) == TRUE)

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
