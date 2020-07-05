# TODO: Kommentare löschen oder Lösung finden, dass nicht f sondern Name des genutzten kernel angezeigt wird
#lst_kernels <- c(rectangular, triangular, epanechnikov, biweight,
#                 triweight, tricube, gaussian, cosine, logistic,
#                 sigmoidFunction, silverman)

#test_that("the argument of a kernel has to be numeric",
#          sapply(lst_kernels, function(f){expect_error(f("0"))})
#          )

#test_that("the output of a kernel has to be numeric",
#          sapply(lst_kernels, function(f){expect_equal("numeric", class(f(0)))})
#)

#test_that("kernel has to be vectorized",
#          sapply(lst_kernels, function(f){expect_equal(3, length(f(0:2)))})
#)

#test_that("is_kernel has to be true for every kernel",
#          sapply(lst_kernels, function(f){expect_equal(TRUE, is_kernel(f)); print(f)})
#)

test_that("the argument of a kernel has to be numeric",
          c(expect_error(rectangular("0")),
          expect_error(triangular("0")),
          expect_error(epanechnikov("0")),
          expect_error(biweight("0")),
          expect_error(triweight("0")),
          expect_error(tricube("0")),
          expect_error(gaussian("0")),
          expect_error(cosine("0")),
          expect_error(logistic("0")),
          expect_error(sigmoidFunction("0")),
          expect_error(silverman("0")))
)

test_that("the output of a kernel has to be numeric",
          c(expect_equal("numeric", class(rectangular(0))),
          expect_equal("numeric", class(triangular(0))),
          expect_equal("numeric", class(epanechnikov(0))),
          expect_equal("numeric", class(biweight(0))),
          expect_equal("numeric", class(triweight(0))),
          expect_equal("numeric", class(tricube(0))),
          expect_equal("numeric", class(gaussian(0))),
          expect_equal("numeric", class(cosine(0))),
          expect_equal("numeric", class(logistic(0))),
          expect_equal("numeric", class(sigmoidFunction(0))),
          expect_equal("numeric", class(silverman(0))))
)

test_that("kernel has to be vectorized",
          c(expect_equal(3, length(rectangular(0:2))),
          expect_equal(3, length(triangular(0:2))),
          expect_equal(3, length(epanechnikov(0:2))),
          expect_equal(3, length(biweight(0:2))),
          expect_equal(3, length(triweight(0:2))),
          expect_equal(3, length(tricube(0:2))),
          expect_equal(3, length(gaussian(0:2))),
          expect_equal(3, length(cosine(0:2))),
          expect_equal(3, length(logistic(0:2))),
          expect_equal(3, length(sigmoidFunction(0:2))),
          expect_equal(3, length(silverman(0:2))))
)

test_that("is_kernel has to be true for every kernel",
          c(expect_equal(TRUE, is_kernel(rectangular)),
          expect_equal(TRUE, is_kernel(triangular)),
          expect_equal(TRUE, is_kernel(epanechnikov)),
          expect_equal(TRUE, is_kernel(biweight)),
          expect_equal(TRUE, is_kernel(triweight)),
          expect_equal(TRUE, is_kernel(tricube)),
          expect_equal(TRUE, is_kernel(gaussian)),
          expect_equal(TRUE, is_kernel(cosine)),
          expect_equal(TRUE, is_kernel(logistic)),
          expect_equal(TRUE, is_kernel(sigmoidFunction)),
          expect_equal(TRUE, is_kernel(silverman)))
)
