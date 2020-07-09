#'Rejection Sampling (I.e., Accept-reject Method)
#'
#'The \code{rejection_sampling} function operator creates a
#'function that draws samples from a difficult probability density function.
#'
#'@param f_den S3 object of the class \link[KDE:Density]{Density}. It's function
#'  \code{f_den$fun} is a normalized probability density function, from which
#'  the desired samples will be drawn.
#'@param g_den S3 object of the class \link[KDE:Density]{Density}. It's function
#'  \code{g_den$fun} is a normalized probability density function that is used
#'  as a helper-function to draw samples from, instead of drawing samples from
#'  \code{f_den}.
#'@param g random number generator function that draws samples from the density
#'  function \code{g_den$fun}.
#'@param M real number, which satisfies the condition that \code{f(x)} is less
#'  or equal than \code{M*g(x)} for all real numbers x.
#'@details the algorithm first draws samples from the more "well-behaved"
#'  density function \code{g_den$fun} and then accepts/rejects these draws,
#'  according to whether they are likely within the proposed denisty function
#'  \code{f_den$fun} or not. It is the most efficient if \code{M} is chosen as small
#'  as possible. That means that the shape of g_den is very close to f_den and
#'  as a result it is more likely that the samples drawn from \code{g_den$fun} get
#'  accepted.
#' @examples
#' custom_den <- function(x) {
#' ret <- 1 + sin(2*pi*x)
#' ret[x < 0 | 1 < x] <- 0
#' ret
#' }
#' f_den <- Density(custom_den, c(0,1))
#' g_den <- Density(dunif)
#'
#'custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)
#'x <- seq(-0.5, 1.5, by=0.01)
#'y <- f_den$fun(x)
#'
#'plot(x, y, type="l", main="Custom density: 1 + sin(2*pi*x)", ylab="density")
#'n <- 65
#'points(custom_sampler(65), rep(1, n), col="red", cex=0.5)
#'
#'@return returning a function with argument \code{n} that draws \code{n}
#'  samples from the density function \code{f_den$fun}.
#'@export

rejection_sampling <- function(f_den, g_den, g, M){
  # M is constant in the real numbers
  stopifnot(is.numeric(M))
  stopifnot("M has length 1" = all.equal(length(M), 1))

  # f_den, g_den should be probability density functions
  stopifnot("f_den has to be a Density obj" = all.equal(validate_Density(f_den), f_den))
  stopifnot("g_den has to be a Density obj" = all.equal(validate_Density(g_den), g_den))
  stopifnot("g has to be a numerical function" = class(g) == "function")
  # relation between f_den,M and g_den that has to be satisfied
  lower <- f_den$support[1]
  upper <- f_den$support[2]

  testing_points <- seq(max(lower, -1e10), min(upper, 1e10), length.out=1e4)

  stopifnot("relation between f_den,M and g_den that has to be satisfied"= all(f_den$fun(testing_points) <= M * g_den$fun(testing_points)))

  sampling(f_den$fun, g_den$fun, g, M)
  }

sampling <- function(f_den, g_den, g, M) {
  force(f_den)
  force(g_den)
  force(g)
  force(M)

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
