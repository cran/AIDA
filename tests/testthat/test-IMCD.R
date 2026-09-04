library(testthat)

test_that("IMCD returns raw estimates when m equals n", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")

  res <- IMCD(obj, m = obj@NObs, cutoff = "raw")

  expect_equal(res$cutoff, "raw")
  expect_true(all(res$final_z == 1))
  expect_true(is.matrix(res$cov_IMCD))
  expect_equal(res$cov_IMCD, int_cov_z(rep(1, obj@NObs), obj))
  expect_equal(res$mean_IMCD_c, int_mean_z(rep(1, obj@NObs), as.matrix(obj@Centers)))
  expect_equal(res$mean_IMCD_r, int_mean_z(rep(1, obj@NObs), as.matrix(obj@Ranges)))
  expect_equal(res$cutoff_value, NA)
  expect_true(all(names(res$robust_dist) == rownames(obj)))
})

test_that("IMCD errors on non-intData input", {
  expect_error(IMCD(1:3), "Argument data is not an object of class intData")
})

test_that("IMCD supports adjbox, F-dist and farness cutoffs when dependencies are available", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")

  # adjbox (robustbase)
  testthat::skip_if_not_installed("robustbase")
  res_adj <- IMCD(obj, m = floor(obj@NObs * 0.75), cutoff = "adjbox")
  expect_equal(res_adj$cutoff, "adjbox")
  expect_true(!is.null(res_adj$cutoff_value))

  # F-dist (CerioliOutlierDetection)
  testthat::skip_if_not_installed("CerioliOutlierDetection")
  res_F <- IMCD(obj, m = floor(obj@NObs * 0.75), cutoff = "F-dist")
  expect_equal(res_F$cutoff, "F-dist")
  expect_true(!is.null(res_F$cutoff_value))

  # farness (cellWise)
  testthat::skip_if_not_installed("cellWise")
  res_far <- IMCD(obj, m = floor(obj@NObs * 0.75), cutoff = "farness")
  expect_equal(res_far$cutoff, "farness")
  expect_true(!is.null(res_far$cutoff_value) || is.na(res_far$cutoff_value))
})

test_that("IMCD supports data with 1 variable", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6, 7, 8),
    U1 = c(3, 4, 5, 6, 7, 8, 9, 10)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "X",
                 LatentParam = list(0.25), LatentCase = "U_id_symmetric",
                 LatentDist = "Unif")

  res <- IMCD(obj, cutoff = "raw")

  expect_true(is.matrix(res$cov_IMCD))
  expect_equal(dim(res$cov_IMCD), c(1, 1))
  expect_length(res$mean_IMCD_c, 1)
  expect_length(res$mean_IMCD_r, 1)
})

test_that("IMCD throws error when m<=p", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  
  expect_error(IMCD(obj, m = obj@NIVar), "Subset size m must satisfy m > p for positive definiteness.")
})

test_that("IMCD defaults m to floor(0.75 * n) when m is NULL", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  res <- IMCD(obj, cutoff = "raw")
  expect_equal(sum(res$final_z), floor(0.75 * obj@NObs))
})

test_that("IMCD uses bigIMCD for large samples and returns correct chi-squared cutoff", {
  set.seed(2027)
  n <- 610
  Data <- data.frame(L1 = rnorm(n), U1 = rnorm(n, 1), L2 = rnorm(n), U2 = rnorm(n, 1))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("A","B"))
  res <- IMCD(obj, m = floor(0.75 * n), cutoff = "chi-squared", cutoff_lvl = 0.5)
  expect_equal(res$cutoff, "chi-squared")

  expect_equal(res$cutoff_value, qchisq(0.5, df = obj@NIVar))
  expect_equal(length(res$robust_dist), obj@NObs)
  expect_true(is.numeric(res$robust_dist))
})
