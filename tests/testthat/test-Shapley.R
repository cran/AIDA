library(testthat)

setup_shapley_data <- function() {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4, 5, 6), U1 = c(3, 4, 5, 6, 7, 8),
    L2 = c(0, 1, 0, -1, -1, 1), U2 = c(2, 3, 3, 2, 3, 2)
  )
  intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
          LatentParam = list(0.25, 0.5), LatentCase = "U_id", LatentDist = "Unif")
}

test_that("int_Shapley returns matrix with observation and variable names", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)

  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)

  expect_true(is.matrix(shapley))
  expect_equal(dim(shapley), c(obj@NObs, obj@NIVar))
  expect_equal(rownames(shapley), rownames(obj))
  expect_equal(colnames(shapley), colnames(obj))
})

test_that("int_Shapley values sum to the Interval-Mahalanobis distance", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)

  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  d2 <- IMah_dist(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)

  expect_equal(as.numeric(rowSums(shapley)), as.numeric(d2), tolerance = 1e-8)
})

test_that("int_Shapley_decomp returns a list of matrices with expected components", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)

  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  d2 <- IMah_dist(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)

  expect_type(decomposed, "list")
  expect_equal(length(decomposed), obj@NObs)
  expect_equal(names(decomposed), rownames(obj))
  for (i in seq_len(obj@NObs)) {
    mat <- decomposed[[i]]
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), c(3, obj@NIVar))
    expect_equal(rownames(mat), c("Centers", "Ranges", "CentersRanges"))
    expect_equal(colnames(mat), colnames(obj))
    expect_equal(sum(mat), sum(shapley[i, ]), tolerance = 1e-8)
    expect_equal(sum(mat), d2[[i]], tolerance = 1e-8)
  }
})

test_that("int_Shapley_interaction returns symmetric matrices with variable names", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)

  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  d2 <- IMah_dist(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)

  expect_type(interaction, "list")
  expect_equal(length(interaction), obj@NObs)
  expect_equal(names(interaction), rownames(obj))
  for (i in seq_len(obj@NObs)) {
    mat <- interaction[[i]]
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), c(obj@NIVar, obj@NIVar))
    expect_equal(rownames(mat), colnames(obj))
    expect_equal(colnames(mat), colnames(obj))
    expect_true(isTRUE(all.equal(mat, t(mat), tolerance = 1e-8)))
    expect_equal(sum(mat), sum(shapley[i, ]), tolerance = 1e-8)
    expect_equal(sum(mat), d2[[i]], tolerance = 1e-8)
  }
})

test_that("Shapley functions error on non-intData input", {
  expect_error(int_Shapley(1:3), "Argument data is not an object of class intData")
  expect_error(int_Shapley_decomp(1:3), "Argument data is not an object of class intData")
  expect_error(int_Shapley_interaction(1:3), "Argument data is not an object of class intData")
})

test_that("int_Shapley works with U_id_symmetric LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4), U1 = c(3, 4, 5, 6),
    L2 = c(0, 1, 0, -1), U2 = c(2, 3, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_true(is.matrix(shapley))
  expect_equal(dim(shapley), c(obj@NObs, obj@NIVar))
  expect_false(any(is.na(shapley)))
  expect_false(any(is.infinite(shapley)))
})

test_that("int_Shapley_decomp works with U_id_symmetric LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3), U1 = c(3, 4, 5),
    L2 = c(0, 1, 0), U2 = c(2, 3, 3)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_type(decomposed, "list")
  expect_equal(length(decomposed), obj@NObs)
  for (i in seq_len(obj@NObs)) {
    mat <- decomposed[[i]]
    expect_equal(dim(mat), c(2, obj@NIVar))
    expect_equal(rownames(mat), c("Centers", "Ranges"))
  }
})

test_that("int_Shapley_interaction works with U_id_symmetric LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3), U1 = c(3, 4, 5),
    L2 = c(0, 1, 0), U2 = c(2, 3, 3)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_type(interaction, "list")
  expect_equal(length(interaction), obj@NObs)
  for (i in seq_len(obj@NObs)) {
    mat <- interaction[[i]]
    expect_true(is.matrix(mat))
    expect_true(isTRUE(all.equal(mat, t(mat), tolerance = 1e-8)))
  }
})

test_that("int_Shapley works with General LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4), U1 = c(3, 4, 5, 6),
    L2 = c(0, 1, 0, -1), U2 = c(2, 3, 3, 2)
  )
  e_UU <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
  psi <- diag(c(0.3, 0.3))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(e_UU, psi), LatentCase = "General", LatentDist = "Normal")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_true(is.matrix(shapley))
  expect_equal(dim(shapley), c(obj@NObs, obj@NIVar))
  expect_false(any(is.na(shapley)))
  expect_false(any(is.infinite(shapley)))
})

test_that("int_Shapley_decomp works with General LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3), U1 = c(3, 4, 5),
    L2 = c(0, 1, 0), U2 = c(2, 3, 3)
  )
  e_UU <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
  psi <- diag(c(0.3, 0.3))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(e_UU, psi), LatentCase = "General", LatentDist = "Normal")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_type(decomposed, "list")
  for (i in seq_len(obj@NObs)) {
    mat <- decomposed[[i]]
    expect_equal(dim(mat), c(3, obj@NIVar))
    expect_equal(rownames(mat), c("Centers", "Ranges", "CentersRanges"))
  }
})

test_that("int_Shapley_interaction works with General LatentCase", {
  Data <- data.frame(
    L1 = c(1, 2, 3), U1 = c(3, 4, 5),
    L2 = c(0, 1, 0), U2 = c(2, 3, 3)
  )
  e_UU <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
  psi <- diag(c(0.3, 0.3))
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(e_UU, psi), LatentCase = "General", LatentDist = "Normal")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_type(interaction, "list")
  for (i in seq_len(obj@NObs)) {
    mat <- interaction[[i]]
    expect_true(isTRUE(all.equal(mat, t(mat), tolerance = 1e-8)))
  }
})

test_that("int_Shapley computes missing parameters with IMCD", {
  obj <- setup_shapley_data()
  
  # Without any parameters
  shapley_full <- int_Shapley(obj)
  expect_true(is.matrix(shapley_full))
  expect_equal(dim(shapley_full), c(obj@NObs, obj@NIVar))
})

test_that("int_Shapley_decomp computes missing parameters with IMCD", {
  obj <- setup_shapley_data()
  
  # Without any parameters
  decomposed_full <- int_Shapley_decomp(obj)
  expect_type(decomposed_full, "list")
  expect_equal(length(decomposed_full), obj@NObs)
})

test_that("int_Shapley_interaction computes missing parameters with IMCD", {
  obj <- setup_shapley_data()
  
  # Without any parameters
  interaction_full <- int_Shapley_interaction(obj)
  expect_type(interaction_full, "list")
  expect_equal(length(interaction_full), obj@NObs)
})

test_that("int_Shapley works with single observation", {
  Data <- data.frame(
    L1 = c(1), U1 = c(3),
    L2 = c(0), U2 = c(2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25, 0.5), LatentCase = "U_id", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_equal(dim(shapley), c(1, 2))
  expect_equal(rownames(shapley), rownames(obj))
})

test_that("int_Shapley works with multiple variables", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 3, 4, 5), U1 = c(3, 4, 5, 1, 2, 3),
    L2 = c(0, 1, 0, 1, 2, 3), U2 = c(2, 3, 3, 3, 4, 5),
    L3 = c(1, 2, 3, 0, 1, 0), U3 = c(2, 3, 4, 2, 3, 3)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y", "Z"),
                 LatentParam = list(0.25, 0.5), LatentCase = "U_id", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_equal(dim(shapley), c(6, 3))
  for (i in seq_len(obj@NObs)) {
    expect_equal(dim(interaction[[i]]), c(3, 3))
  }
})

test_that("int_Shapley_decomp components sum correctly for U_id", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  for (i in seq_len(obj@NObs)) {
    mat <- decomposed[[i]]
    expected_shapley <- colSums(mat)
    actual_shapley <- shapley[i, ]
    expect_equal(as.numeric(expected_shapley), as.numeric(actual_shapley), tolerance = 1e-8)
  }
})

# Verify decomposition components for U_id_symmetric sum correctly
test_that("int_Shapley_decomp components sum correctly for U_id_symmetric", {
  Data <- data.frame(
    L1 = c(1, 2, 3, 4), U1 = c(3, 4, 5, 6),
    L2 = c(0, 1, 0, -1), U2 = c(2, 3, 3, 2)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25), LatentCase = "U_id_symmetric", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  for (i in seq_len(obj@NObs)) {
    mat <- decomposed[[i]]
    expected_shapley <- colSums(mat)
    actual_shapley <- shapley[i, ]
    expect_equal(as.numeric(expected_shapley), as.numeric(actual_shapley), tolerance = 1e-8)
  }
})

# Verify interaction matrix diagonal properties
test_that("int_Shapley_interaction matrix diagonal and sum properties", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  cov_inv <- safe_solve_cov(cov)
  
  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  for (i in seq_len(obj@NObs)) {
    mat <- interaction[[i]]
    diag_values <- diag(mat)

    # Compute expected diagonal elements
    c_0 <- C[i,]-mean_c
    r_0 <- R[i,]-mean_r
    delta <- obj@LatentParam[[1]]
    U_mean <- obj@LatentParam[[2]]
    inter_shapley <- 2*tcrossprod(c_0)*cov_inv + 2*delta*tcrossprod(r_0)*cov_inv + 
                    + U_mean*tcrossprod(c_0,r_0)*cov_inv + U_mean*tcrossprod(r_0,c_0)*cov_inv

    expected_diag <- diag(inter_shapley - diag(shapley[i, ]))
    expect_equal(as.numeric(diag_values), as.numeric(expected_diag), tolerance = 1e-8)
  }
})

# Test variable names preservation through all functions
test_that("Variable and row names preserved through Shapley functions", {
  obj <- setup_shapley_data()
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  decomposed <- int_Shapley_decomp(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  interaction <- int_Shapley_interaction(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  # Check shapley names
  expect_equal(colnames(shapley), colnames(obj))
  expect_equal(rownames(shapley), rownames(obj))
  
  # Check decomposed names
  for (i in seq_len(obj@NObs)) {
    expect_equal(colnames(decomposed[[i]]), colnames(obj))
  }
  
  # Check interaction names
  for (i in seq_len(obj@NObs)) {
    expect_equal(rownames(interaction[[i]]), colnames(obj))
    expect_equal(colnames(interaction[[i]]), colnames(obj))
  }
})

# Test numerical stability with extreme values
test_that("int_Shapley handles different scales of data", {
  Data <- data.frame(
    L1 = c(100, 200, 300), U1 = c(300, 400, 500),
    L2 = c(0.001, 0.002, 0.003), U2 = c(0.01, 0.02, 0.03)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25, 0.5), LatentCase = "U_id", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_false(any(is.na(shapley)))
  expect_false(any(is.infinite(shapley)))
})

# Test with negative values
test_that("int_Shapley handles negative interval values", {
  Data <- data.frame(
    L1 = c(-5, -2, 0), U1 = c(-1, 2, 5),
    L2 = c(-3, -1, 1), U2 = c(-1, 1, 3)
  )
  obj <- intData(Data, Seq = "LbUb_VarbyVar", VarNames = c("X", "Y"),
                 LatentParam = list(0.25, 0.5), LatentCase = "U_id", LatentDist = "Unif")
  
  C <- as.matrix(obj@Centers)
  R <- as.matrix(obj@Ranges)
  mean_c <- colMeans(C)
  mean_r <- colMeans(R)
  cov <- int_cov(obj)
  
  shapley <- int_Shapley(obj, mean_c = mean_c, mean_r = mean_r, cov = cov)
  
  expect_equal(dim(shapley), c(3, 2))
  expect_false(any(is.na(shapley)))
})
