library(testthat)

test_that("int_outliers with chi-squared cutoff returns expected structure", {
  set.seed(1)
  robust_dist <- abs(rnorm(6))
  names(robust_dist) <- paste0("obs", seq_along(robust_dist))

  res <- int_outliers(robust_dist, cutoff = "chi-squared", p = 2)

  expect_type(res, "list")
  expect_true("is_outlier" %in% names(res))
  expect_equal(length(res$is_outlier), length(robust_dist))
  expect_equal(res$cutoff, "chi-squared")
  expect_equal(res$cutoff_value, qchisq(0.975, df = 2))
})

test_that("int_outliers errors when p is missing for chi-squared and F-dist", {
  robust_dist <- abs(rnorm(5))
  expect_error(int_outliers(robust_dist, cutoff = "chi-squared"), "must provide the number of variables p")
  expect_error(int_outliers(robust_dist, cutoff = "F-dist"), "must provide the number of variables p")
})

test_that("int_outliers adjbox returns fence when robustbase installed", {
  testthat::skip_if_not_installed("robustbase")
  robust_dist <- abs(rnorm(10))
  names(robust_dist) <- paste0("r", seq_along(robust_dist))

  res <- int_outliers(robust_dist, cutoff = "adjbox")
  expect_equal(res$cutoff, "adjbox")
  expect_true(is.numeric(res$cutoff_value))
  expect_true(length(res$cutoff_value) == 2)
  expect_equal(length(res$is_outlier), length(robust_dist))
})

test_that("int_outliers F-dist requires z and p and returns expected fields when provided", {
  testthat::skip_if_not_installed("CerioliOutlierDetection")
  robust_dist <- abs(rnorm(8))
  names(robust_dist) <- paste0("o", seq_along(robust_dist))
  z <- rep(c(1,0), length.out = length(robust_dist))

  # provide minimal valid inputs for F-dist
  expect_error(int_outliers(robust_dist, cutoff = "F-dist", p = NULL, z = NULL))
  # When z and p provided, function should run (may rely on qbeta/qf internals)
  res <- int_outliers(robust_dist, cutoff = "F-dist", p = 2, z = z)
  expect_true(is.list(res))
  expect_true("is_outlier" %in% names(res))
})

test_that("int_outliers farness uses cellWise when available", {
  testthat::skip_if_not_installed("cellWise")
  robust_dist <- abs(rnorm(12))
  names(robust_dist) <- paste0("a", seq_along(robust_dist))

  res <- int_outliers(robust_dist, cutoff = "farness", cutoff_lvl = 0.9)
  expect_equal(res$cutoff, "farness")
  expect_true(!is.null(res$farness_probs))
  expect_equal(length(res$is_outlier), length(robust_dist))
})

test_that("int_outliers flags a clear extreme as outlier (chi-squared)", {
  robust_dist <- c(a = 0.1, b = 0.2, extreme = 100)
  res <- int_outliers(robust_dist, cutoff = "chi-squared", p = 2)

  expect_true(res$is_outlier["extreme"])
  expect_equal(res$outliers_names, "extreme")
  expect_equal(names(res$is_outlier), names(robust_dist))
})

test_that("int_outliers adjbox respects cutoff_lvl when robustbase available", {
  testthat::skip_if_not_installed("robustbase")
  robust_dist <- c(1, 2, 3, 10, 20)
  names(robust_dist) <- paste0("x", seq_along(robust_dist))

  res <- int_outliers(robust_dist, cutoff = "adjbox", cutoff_lvl = 2)
  expect_equal(res$cutoff, "adjbox")
  expected_fence <- robustbase::adjboxStats(robust_dist, coef = 2, doScale = FALSE)$fence
  expect_equal(res$cutoff_value, expected_fence)
  expect_equal(length(res$is_outlier), length(robust_dist))
})
