test_that("kernels have to work on validate_Kernel",
          c(expect_equal(rectangular, validate_Kernel(rectangular)),
            expect_equal(triangular, validate_Kernel(triangular)),
            expect_equal(epanechnikov, validate_Kernel(epanechnikov)),
            expect_equal(biweight, validate_Kernel(biweight)),
            expect_equal(triweight, validate_Kernel(triweight)),
            expect_equal(tricube, validate_Kernel(tricube)),
            expect_equal(gaussian, validate_Kernel(gaussian)),
            expect_equal(cosine, validate_Kernel(cosine)),
            expect_equal(logistic, validate_Kernel(logistic)),
            expect_equal(sigmoid, validate_Kernel(sigmoid)),
            expect_equal(silverman, validate_Kernel(silverman)))
)

test_that("the argument of a kernel has to be numeric",
          c(expect_error(rectangular$fun("0")),
          expect_error(triangular$fun("0")),
          expect_error(epanechnikov$fun("0")),
          expect_error(biweight$fun("0")),
          expect_error(triweight$fun("0")),
          expect_error(tricube$fun("0")),
          expect_error(gaussian$fun("0")),
          expect_error(cosine$fun("0")),
          expect_error(logistic$fun("0")),
          expect_error(sigmoid$fun("0")),
          expect_error(silverman$fun("0")))
)

# all built-in kernels have a local max in zero
test_that("the built-in functions are centered in zero",{
          c(expect_equal("numeric", class(rectangular$fun(0))),
          expect_equal("numeric", class(triangular$fun(0))),
          expect_equal("numeric", class(epanechnikov$fun(0))),
          expect_equal("numeric", class(biweight$fun(0))),
          expect_equal("numeric", class(triweight$fun(0))),
          expect_equal("numeric", class(tricube$fun(0))),
          expect_equal("numeric", class(gaussian$fun(0))),
          expect_equal("numeric", class(cosine$fun(0))),
          expect_equal("numeric", class(logistic$fun(0))),
          expect_equal("numeric", class(sigmoid$fun(0))),
          expect_equal("numeric", class(silverman$fun(0))))
          x <- seq(-1, 1, length.out=1000)
          x <- x[x!=0]
          expect_true(all(rectangular$fun(x) <= rectangular$fun(0)))
          expect_true(all(triangular$fun(x) <= triangular$fun(0)))
          expect_true(all(epanechnikov$fun(x) <= epanechnikov$fun(0)))
          expect_true(all(biweight$fun(x) <= biweight$fun(0)))
          expect_true(all(triweight$fun(x) <= triweight$fun(0)))
          expect_true(all(tricube$fun(x) <= tricube$fun(0)))
          expect_true(all(gaussian$fun(x) <= gaussian$fun(0)))
          expect_true(all(cosine$fun(x) <= cosine$fun(0)))
          expect_true(all(logistic$fun(x) <= logistic$fun(0)))
          expect_true(all(sigmoid$fun(x) <= sigmoid$fun(0)))
          expect_true(all(silverman$fun(x) <= silverman$fun(0)))
  }
)


test_that("kernel has to be vectorized",
          c(expect_equal(3, length(rectangular$fun(0:2))),
          expect_equal(3, length(triangular$fun(0:2))),
          expect_equal(3, length(epanechnikov$fun(0:2))),
          expect_equal(3, length(biweight$fun(0:2))),
          expect_equal(3, length(triweight$fun(0:2))),
          expect_equal(3, length(tricube$fun(0:2))),
          expect_equal(3, length(gaussian$fun(0:2))),
          expect_equal(3, length(cosine$fun(0:2))),
          expect_equal(3, length(logistic$fun(0:2))),
          expect_equal(3, length(sigmoid$fun(0:2))),
          expect_equal(3, length(silverman$fun(0:2))))
)

