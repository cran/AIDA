library(testthat)

test_that("Mallows_dist computes squared distances for U_id_symmetric", {
  Data <- data.frame(L1 = c(1, 2), U1 = c(3, 4))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V", LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  result <- Mallows_dist(obj, mean_c = 2.5, mean_r = 2)

  expect_named(result, c("1", "2"))
  expect_equal(result, c('1'=0.25, '2'=0.25))
})

test_that("Mallows_dist computes the U_id cross term correctly", {
  Data <- data.frame(L1 = c(1, 3), U1 = c(4, 7))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V", LatentParam = list(0.5, 0.2), LatentCase = "U_id", LatentDist = "Triang")
  result <- Mallows_dist(obj, mean_c = 3.75, mean_r = 3.5)

  expect_equal(result, c('1'=1.8125, '2'=1.8125))
})

test_that("Mallows_dist computes the General case using automatic means", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  Param <- list(matrix(c(0.5,0.25,0.12,0.5),2,2), matrix(c(1,0,0,1),2,2))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X","Y"),
                 LatentParam = Param,
                 LatentCase = "General", LatentDist = "KDE")

  result <- Mallows_dist(obj)

  # Compute expected result
  mean_c <- colMeans(obj@Centers)
  mean_r <- colMeans(obj@Ranges)
  C_0 <- t(obj@Centers)-mean_c
  R_0 <- t(obj@Ranges)-mean_r

  delta <- 1/4*diag(diag(Param[[1]]))
  psi <- Param[[2]]

  expected <- diag(crossprod(C_0)) + 
                + diag(crossprod(R_0,delta%*%R_0)) +
                + diag(crossprod(C_0,psi%*%R_0))

  names(expected) <- rownames(obj)

  expect_named(result, rownames(obj))
  expect_equal(result, expected)
})

test_that("IMah_dist handles z-based subsetting and robust covariance fallback", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, 1, 0, 1), U2 = c(2, 3, 2, 3, 2, 3)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")

  z <- c(1, 0, 1, 1, 0, 1)
  result <- IMah_dist(obj, z = z)

  # Computed expected result
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  param <- obj@LatentParam
  case <- obj@LatentCase

  mean_c <- int_mean_z(z, C)
  mean_r <- int_mean_z(z, R)
  cov <- int_cov_z(z,obj)
  cov_inv <- safe_solve_cov(cov)
  delta <- param[[1]]
  expected <- mahalanobis(C, mean_c, cov_inv, inverted = TRUE) + 
              + delta*mahalanobis(R, mean_r, cov_inv, inverted = TRUE)
  names(expected) <- rownames(obj)

  expect_named(result, rownames(obj))
  expect_equal(result, expected)
})

test_that("IMah_dist computes the U_id case with explicit covariance", {
  Data <- data.frame(L1 = c(1, 3), U1 = c(4, 7))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V",
                 LatentParam = list(0.5, 0.2), LatentCase = "U_id", LatentDist = "Triang")
  cov_mat <- diag(1, 1)
  mean_c <- 3.75
  mean_r <- 3.5

  result <- IMah_dist(obj, mean_c = mean_c, mean_r = mean_r, cov = cov_mat)

  # Compute expected result
  delta <- 0.5
  U_mean <- 0.2
  C_0 <- t(obj@Centers)-mean_c
  R_0 <- t(obj@Ranges)-mean_r
  cov_inv <- safe_solve_cov(cov_mat)
  expected <- mahalanobis(obj@Centers, mean_c, cov_inv, inverted = TRUE) + 
               + delta*mahalanobis(obj@Ranges, mean_r, cov_inv, inverted = TRUE) +
               + U_mean*diag(crossprod(C_0,cov_inv%*%R_0))

  expect_equal(result, expected)
})

test_that("IMah_dist_pairs computes U_id pairwise distances correctly", {
  Data <- data.frame(L1 = c(1, 3), U1 = c(4, 7))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V",
                 LatentParam = list(0.5, 0.2), LatentCase = "U_id", LatentDist = "Triang")
  cov_mat <- diag(1, 1)

  result <- IMah_dist_pairs(obj, cov = cov_mat)
  expected <- matrix(c(0, 7.25, 7.25, 0), nrow = 2)
  colnames(expected) <- rownames(expected) <- c("1", "2")

  expect_equal(result, expected)
})

test_that("IMah_dist_pairs computes the General case with explicit covariance", {
  Data <- data.frame(L1 = c(1, 3), U1 = c(2, 5))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V",
                 LatentParam = list(matrix(0.5,1,1), matrix(2,1,1)),
                 LatentCase = "General", LatentDist = "KDE")
  cov_mat <- diag(1, 1)

  result <- IMah_dist_pairs(obj, cov = cov_mat)
  expected <- matrix(c(0, 11.375, 11.375, 0), nrow = 2)
  colnames(expected) <- rownames(expected) <- c("1", "2")

  expect_equal(result, expected)
})

test_that("IMah_dist uses explicit mean and covariance", {
  Data <- data.frame(L1 = c(1, 2), U1 = c(3, 5))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V", LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  cov_mat <- diag(1, 1)
  mean_c <- 2
  mean_r <- 2

  result <- IMah_dist(obj, mean_c = mean_c, mean_r = mean_r, cov = cov_mat)
  expected <- mahalanobis(obj@Centers, mean_c, solve(cov_mat), inverted = TRUE) + 
                + 0.25*mahalanobis(obj@Ranges, mean_r, solve(cov_mat), inverted = TRUE)
  names(expected) <- c("1", "2")
  expect_equal(result, expected)
})

test_that("IMah_dist_pairs computes symmetric pairwise distances", {
  Data <- data.frame(L1 = c(1, 3), U1 = c(3, 6))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = "V", LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  cov_mat <- diag(1, 1)

  result <- IMah_dist_pairs(obj, cov = cov_mat)
  expected_offdiag <- 6.5
  expected <- matrix(c(0, expected_offdiag, expected_offdiag, 0), nrow = 2)
  colnames(expected) <- rownames(expected) <- c("1", "2")

  expect_equal(result, expected)
  expect_equal(result, t(result))
})

test_that("distance functions error on non-intData input", {
  expect_error(Mallows_dist(1:3), "Argument data is not an object of class intData")
  expect_error(IMah_dist(1:3), "Argument data is not an object of class intData")
  expect_error(IMah_dist_pairs(1:3), "Argument data is not an object of class intData")
})
