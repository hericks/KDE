test_that("a kernel should be recognized as a kernel",
          expect_true(is_kernel(gaussian))
)

transformed_kernel <- kernelTransform(gaussian, 2, 2)
# integrate does not integrate perfectly -> deviation from 1 < integration error
transformed_kernel2 <- kernelTransform(rectangular, 0, 0.01)
# function with a. e. zero will crash
a_e_zero <- kernelTransform(rectangular, 0, 0.001)
test_that("a transformed kernel is still a kernel",
          c(expect_true(is_kernel(transformed_kernel)),
          expect_true(is_kernel(transformed_kernel2)),
          expect_true(is_kernel(a_e_zero)))
)

div_integral_fun <- function(x) 1
finite_integral_fun <- function(x) (x <= 10)
non_integrable_fun <- function(x) 1/x
test_that("the integral of a kernel over the real numbers equals one",
          c(expect_false(is_kernel(div_integral_fun)),
            expect_false(is_kernel(finite_integral_fun)),
            expect_false(is_kernel(non_integrable_fun))
            )
)

non_numeric_fun <- function(x){"x"}
test_that("kernel is a numeric function",
          c(expect_false(is_kernel("hallo")),
          expect_false(is_kernel(non_numeric_fun)))
          )
