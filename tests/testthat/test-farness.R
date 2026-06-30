library(testthat)

test_that("farness returns probs vector without cutoff_value", {
  testthat::skip_if_not_installed("cellWise")
  set.seed(42)
  dist <- abs(rnorm(10))

  probs <- farness(dist)

  expect_type(probs, "double")
  expect_equal(length(probs), length(dist))
  expect_true(all(probs >= 0 & probs <= 1))
})

test_that("farness with cutoff_value returns list with farness_probs and cutoff_value", {
  testthat::skip_if_not_installed("cellWise")
  set.seed(42)
  dist <- abs(rnorm(8))

  res <- farness(dist, cutoff_value = 0.9)

  expect_type(res, "list")
  expect_true("farness_probs" %in% names(res))
  expect_true("cutoff_value" %in% names(res))
  expect_equal(length(res$farness_probs), length(dist))
  expect_true(is.numeric(res$cutoff_value))
})

test_that("farness errors when cutoff_value is outside [0, 1]", {
  testthat::skip_if_not_installed("cellWise")
  dist <- abs(rnorm(100))

  expect_error(farness(dist, cutoff_value = -0.1), "cutoff_value must be between 0 and 1")
  expect_error(farness(dist, cutoff_value = 1.5), "cutoff_value must be between 0 and 1")
})

test_that("farness handles named distances", {
  testthat::skip_if_not_installed("cellWise")
  dist <- abs(rnorm(100))
  names(dist) <- rep(c("a","b","c","d","e"), 20)
  
  probs <- farness(dist)

  expect_equal(names(probs), names(dist))
  expect_equal(length(probs), length(dist))
})

test_that("farness filters out very small distances", {
  testthat::skip_if_not_installed("cellWise")
  dist <- c(rep(c(1e-12, 1e-11), 20), abs(rnorm(100)))

  probs <- farness(dist)

  expect_equal(length(probs), length(dist))
  expect_true(all(!is.na(probs)))
})

test_that("farness with cutoff_value=0.5 returns valid cutoff in original scale", {
  testthat::skip_if_not_installed("cellWise")
  set.seed(123)
  dist <- abs(rnorm(100))

  res <- farness(dist, cutoff_value = 0.5)

  expect_true(is.numeric(res$cutoff_value))
  expect_true(res$cutoff_value > 0)
  expect_true(length(res$cutoff_value) == 1)
})

test_that("farness with large distances returns valid probs", {
  testthat::skip_if_not_installed("cellWise")
  dist <- c(1000, 10000, 100000, 1000000, rnorm(100))

  res <- farness(dist, cutoff_value = 0.95)

  expect_equal(length(res$farness_probs), length(dist))
  expect_true(all(res$farness_probs >= 0 & res$farness_probs <= 1))
  expect_true(!is.na(res$cutoff_value))
})

