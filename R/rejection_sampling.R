#' Rejection Sampling
#'
#' The \code{rejection_sampling} function operator creates a function that draws
#' samples from a given probability density function.
#'
#' @param f_den S3 object of class \code{\link{Density}}; the probability density
#'   to construct the sampler for.
#' @param g_den S3 object of class \code{Density}; the probability density for
#'   the sampler given in \code{g}.
#' @param g \code{R} function with single numeric argument; the random number
#'   generator that draws samples from the density function \code{g_den}.
#' @param M strictly positive numeric scalar; satisfies \code{f(x) <= M*g(x)}
#'   for all numeric scalar inputs \code{x}.
#'
#' @details Rejection sampling uses \code{g} to draw samples and accepts/rejects
#'   these samples according to the densities \code{f_gen} and \code{g_den},
#'   such that the resulting samples are \code{f_den}-distributed. Many rejected
#'   samples result in longer runtimes. To prevent this \code{M} should be
#'   chosen as small as possible, satisfying \code{f_den$fun(x) <=
#'   M*g_den$fun(x)} for all numeric scalar inputs \code{x}.
#'
#' @return A function taking a single numeric scalar argument \code{n},
#'   returning \code{n} \code{f-den} distributed random numbers.
#'
#' @examples
#' custom_den <- function(x) {
#'   ret <- 1 + sin(2*pi*x)
#'   ret[x < 0 | 1 < x] <- 0
#'   ret
#' }
#'
#' f_den <- Density(custom_den, c(0,1))
#' g_den <- Density(dunif)
#'
#' custom_sampler <- rejection_sampling(f_den, g_den, runif, 2)
#' x <- seq(-0.5, 1.5, by=0.01)
#' y <- f_den$fun(x)
#'
#' plot(x, y, type="l", main="Custom density: 1 + sin(2*pi*x)", ylab="density")
#' n <- 65
#' points(custom_sampler(65), rep(0, n), col="red", pch=".", cex=0.5)
#' legend("topright",
#'        legend=c("f_den", "samples"),
#'        col=c("black", "red", "blue"),
#'        pch=c("-", "."))
#'
#' @seealso \code{\link{Density}} for more information about densities.
#'
#' @export
rejection_sampling <- function(f_den, g_den, g, M){
  # M is constant in the real numbers
  stopifnot(is.numeric(M))
  stopifnot("M has length 1" = all.equal(length(M), 1))

  # f_den, g_den should be probability density functions
  stopifnot("f_den has to be a Density obj" = all.equal(validate_Density(f_den), f_den))
  stopifnot("g_den has to be a Density obj" = all.equal(validate_Density(g_den), g_den))
  stopifnot("g has to be a numerical function" = class(g) == "function")
  # g has to be a sampler function
  stopifnot("g has to return the correct number of samples" = length(g(1)) == 1)
  stopifnot("g has to return the correct number of samples" = length(g(10)) == 10)
  stopifnot("g has to return the correct number of samples" = length(g(100)) == 100)
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
