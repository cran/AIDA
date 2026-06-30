test_that("int_mean_z works", {
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), ncol = p)
  z <- c(rep(1, n))
  
  expect_equal(int_mean_z(z, X), colMeans(X))
})

test_that("int_cov computes U_id_symmetric from explicit covariance matrices", {
  sigma_cc <- matrix(c(1, 0, 0, 2), nrow = 2)
  sigma_rr <- matrix(c(2, 0, 0, 3), nrow = 2)
  result <- int_cov(sigma_cc = sigma_cc, sigma_rr = sigma_rr, LatentParam = list(0.25), LatentCase = "U_id_symmetric")
  expect_equal(result, sigma_cc + 0.25 * sigma_rr)
})

test_that("int_cov computes U_id from explicit matrices and sigma_cr", {
  sigma_cc <- matrix(c(1, 1, 1, 2), nrow = 2)
  sigma_rr <- matrix(c(2, 0, 0, 1), nrow = 2)
  sigma_cr <- matrix(c(0.5, -0.5, 0.2, 0.3), nrow = 2)
  latent <- list(0.5, 0.2)
  result <- int_cov(sigma_cc = sigma_cc, sigma_rr = sigma_rr, sigma_cr = sigma_cr, LatentParam = latent, LatentCase = "U_id")
  expected <- sigma_cc + 0.5 * sigma_rr + 0.1 * (sigma_cr + t(sigma_cr))
  expect_equal(result, expected)
})

test_that("int_cov computes General case from explicit matrices", {
  sigma_cc <- matrix(c(1, 2, 2, 3), nrow = 2)
  sigma_rr <- matrix(c(1, 0, 0, 1), nrow = 2)
  sigma_cr <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  e_UU <- matrix(c(1, 0.1, 0.1, 2), nrow = 2)
  psi <- diag(c(0.1, -0.2))
  latent <- list(e_UU, psi)
  result <- int_cov(sigma_cc = sigma_cc, sigma_rr = sigma_rr, sigma_cr = sigma_cr, LatentParam = latent, LatentCase = "General")
  expected <- sigma_cc + 1/4 * (e_UU * sigma_rr) + 1/2 * sigma_cr %*% psi + 1/2 * psi %*% t(sigma_cr)
  expect_equal(result, expected)
})

test_that("int_cov uses intData object and returns named covariance matrix", {
  Data <- data.frame(L1 = c(1, 2), U1 = c(3, 4), L2 = c(0, 1), U2 = c(2, 3))
  VarNames <- c("A", "B")
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = VarNames, LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  result <- int_cov(data = obj)
  expected <- cov(as.matrix(obj@Centers)) + 0.25 * cov(as.matrix(obj@Ranges))
  colnames(expected) <- rownames(expected) <- VarNames
  expect_equal(result, expected)
  expect_equal(colnames(result), VarNames)
  expect_equal(rownames(result), VarNames)
})

test_that("int_cov errors when required arguments are missing", {
  expect_error(int_cov(sigma_cc = matrix(1,1,1), sigma_rr = matrix(1,1,1), LatentParam = list(0.1), LatentCase = "U_id"), "sigma_cr is missing")
  expect_error(int_cov(data = 1:3), "Argument data is not an object of class intData")
})

test_that("int_cov_z computes sample covariance for U_id_symmetric", {
  C <- matrix(c(1, 2, 3, 4), nrow = 2)
  R <- matrix(c(1, 1, 2, 2), nrow = 2)
  data_obj <- intData(cbind(C, C + R), Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  z <- c(1, 0)
  result <- int_cov_z(z, data_obj)
  m <- sum(z)
  C_t <- as.matrix(t(data_obj@Centers))
  R_t <- as.matrix(t(data_obj@Ranges))
  zC <- as.matrix(z * data_obj@Centers)
  zR <- as.matrix(z * data_obj@Ranges)
  sum_c <- C_t %*% zC
  sum_r <- R_t %*% zR
  C_z <- as.matrix(C_t %*% (z %*% t(z)))
  R_z <- as.matrix(R_t %*% (z %*% t(z)))
  expected <- 1/m * sum_c - 1/m^2 * C_z %*% as.matrix(data_obj@Centers) + 0.25/m * sum_r - 0.25/m^2 * R_z %*% as.matrix(data_obj@Ranges)
  colnames(expected) <- rownames(expected) <- c("X", "Y")
  expect_equal(result, expected)
})

test_that("int_cov_z errors for non-intData inputs", {
  expect_error(int_cov_z(c(1, 0), data.frame(x = 1:2)), "Argument data is not an object of class intData")
})

test_that("int_cov_z computes U_id case correctly", {
  # construct small intData with Seq = CenRng_VarbyVar (C1,R1,C2,R2)
  Data <- matrix(c(
    1, 0.5, 2, 0.2,
    3, 0.3, 4, 0.1
  ), nrow = 2, byrow = TRUE)
  obj <- intData(Data, Seq = "CenRng_VarbyVar", VarNames = c("X", "Y"), LatentParam = list(0.25, 0.6), LatentCase = "U_id", LatentDist = "Triang")
  z <- c(1, 1)
  res <- int_cov_z(z, obj)

  # manual computation following the implementation
  param <- obj@LatentParam
  delta <- param[[1]]
  U_mean <- param[[2]]
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  m <- sum(z)
  C_t <- t(C); R_t <- t(R)
  zC <- z * C; zR <- z * R
  sum_c <- C_t %*% zC; sum_r <- R_t %*% zR
  sum_cr <- C_t %*% zR; sum_rc <- R_t %*% zC
  zz_t <- z %*% t(z)
  C_z <- C_t %*% zz_t; R_z <- R_t %*% zz_t

  expected <- 1/m * sum_c - (1/m^2) * C_z %*% C + delta/m * sum_r - (delta/m^2) * R_z %*% R +
    U_mean/(2*m) * sum_cr - U_mean/(2*m^2) * C_z %*% R + U_mean/(2*m) * sum_rc - U_mean/(2*m^2) * R_z %*% C
  colnames(expected) <- rownames(expected) <- c("X", "Y")
  expect_equal(res, expected)
})

test_that("int_cov_z computes General case correctly", {
  Data <- matrix(c(
    1, 0.2, 2, 0.1,
    4, 0.3, 5, 0.4
  ), nrow = 2, byrow = TRUE)
  e_UU <- matrix(c(1, 0.1, 0.1, 2), nrow = 2)
  psi <- diag(c(0.5, -0.2))
  obj <- intData(Data, Seq = "CenRng_VarbyVar", VarNames = c("A", "B"), LatentParam = list(e_UU, psi), LatentCase = "General", LatentDist = "KDE")
  z <- c(1, 0)
  res <- int_cov_z(z, obj)

  # manual computation
  param <- obj@LatentParam
  e_UU_p <- param[[1]]; psi_p <- param[[2]]
  C <- as.matrix(obj@Centers); R <- as.matrix(obj@Ranges); m <- sum(z)
  C_t <- t(C); R_t <- t(R)
  zC <- z * C; zR <- z * R
  sum_c <- C_t %*% zC; sum_r <- R_t %*% zR
  sum_cr <- C_t %*% zR; sum_rc <- R_t %*% zC
  zz_t <- z %*% t(z)
  C_z <- C_t %*% zz_t; R_z <- R_t %*% zz_t

  expected <- 1/m * sum_c - (1/m^2) * C_z %*% C + 1/(4*m) * e_UU_p * sum_r - 1/(4*m^2) * e_UU_p * R_z %*% R +
    1/(2*m) * sum_cr %*% psi_p - 1/(2*m^2) * C_z %*% R %*% psi_p + 1/(2*m) * psi_p %*% sum_rc - 1/(2*m^2) * psi_p %*% R_z %*% C
  colnames(expected) <- rownames(expected) <- c("A", "B")
  expect_equal(res, expected)
})

test_that("safe_solve_cov handles errors, invertible and singular matrices", {
  # non-matrix
  expect_error(safe_solve_cov(1:3), "must be a matrix")
  # non-square
  expect_error(safe_solve_cov(matrix(1:6, nrow = 2)), "must be square")

  # invertible
  cov <- diag(c(2, 3))
  inv_expected <- solve(cov)
  inv_res <- safe_solve_cov(cov, verbose = FALSE)
  expect_equal(inv_res, inv_expected)

  # singular -> fallback to MASS::ginv
  sing <- matrix(c(1, 2, 1, 2), nrow = 2)
  expect_warning(inv_sing <- safe_solve_cov(sing, verbose = TRUE), "Moore-Penrose")
  expect_equal(inv_sing, MASS::ginv(sing))
})
