test_that("integrant should be a function",{
          not_fun <- "x"
          expect_error(integrate_primitive(not_fun, -1, 1))
          a_fun <- function(x) x+1
          expect_error(integrate_primitive(a_fun, -1, 1), NA)
})

test_that("support should make sense",{
  fun <- rectangular$fun
  expect_error(integrate_primitive(fun, rectangular$support[1], rectangular$support[2]), NA)
  expect_error(integrate_primitive(fun, 1, -1))
  expect_error(integrate_primitive(fun, "-1", 1))
  expect_error(integrate_primitive(fun, -1, "-1"))
  expect_error(integrate_primitive(fun, c(-1,0), 1))
  expect_error(integrate_primitive(fun, -1, c(1,2144)))
})

test_that("support should be finite",{
  fun <- rectangular$fun
  expect_error(integrate_primitive(fun, -Inf, -1))
  expect_error(integrate_primitive(fun, -1, Inf))
})

test_that("subdivisions should be numeric value", {
  fun <- rectangular$fun
  expect_error(integrate_primitive(fun, -1, 1, subdivisions=1000L),NA)
  expect_error(integrate_primitive(fun, -1, 1, subdivisions=1000.2342),NA)
  expect_error(integrate_primitive(fun, -1, 1, subdivisions="1000"))
  expect_error(integrate_primitive(fun, -1, 1, subdivisions=c(100, 10)))
})

test_that("non integrable functions will not converge with more subdivisions",{
  non_integrable_fun <- function(x) 1/x
  lower <- 0
  upper <- 1
  rel_error <- c()
  for(i in 1:6){
    subdivisions <- 10^i
    rel_error <- c(rel_error,integrate_primitive(non_integrable_fun, 0, 1, subdivisions, check=TRUE)$relError)
  }
  # rel_error does get smaller using more subdivision steps, but it will not get sufficiently small
  expect_error(stopifnot(res[length(rel_error)] < 0.01))

          })

test_that("integrating integrable functions will converge",{
  convergence <- function(x){
    for(i in seq_along(x)){
      for(j in (i:length(x))){
        if(x[i] < x[j] && i != j) return(FALSE)
      }
    }
    return(TRUE)
  }

  fun <- function(x) x^2 + 2
  lower <- 0
  upper <- 1
  res <- c()
  for(i in 1:5){
    subdivisions <- 10^i
    res <- c(res,integrate_primitive(fun, 0, 1, subdivisions, check=TRUE)$value)

  }
  true_value <- 2+1/3
  expect_true(convergence(abs(res - true_value)))
})

test_that("rel_error has to be boolean value",{
  fun <- function(x) 1/x
  lower <- 0
  upper <- 1
  subdivisions <- 1000L
  res <- c()
  expect_error(integrate_primitive(fun, lower, upper, subdivisions, check=TRUE), NA)
  expect_error(integrate_primitive(fun, lower, upper, subdivisions, check=c(TRUE, FALSE)))
  expect_error(integrate_primitive(fun, lower, upper, subdivisions, check="TRUE"))
  expect_equal(integrate_primitive(fun, lower, upper, subdivisions, check=FALSE)$relError, NULL)
})

