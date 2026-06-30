library(testthat)

test_that("get_latent_var computes expected U and handles bounds", {
  macro <- data.frame(L1 = 1, U1 = 3, L2 = 0, U2 = 2)
  rownames(macro) <- "G1"
  micro <- matrix(c(2,1,
                    3,2,
                    1,0), ncol = 2, byrow = TRUE)
  agrby <- rep("G1", nrow(micro))

  U <- get_latent_var(microdata = micro,
                      macrodata = macro,
                      agrby = agrby,
                      agrlevels = rownames(macro),
                      Seq = "LbUb_VarbyVar")

  centers <- c((macro[1,1] + macro[1,2]) / 2, (macro[1,3] + macro[1,4]) / 2)
  ranges <- c(macro[1,2] - macro[1,1], macro[1,4] - macro[1,3])

  expected <- matrix(NA_real_, nrow = nrow(micro), ncol = ncol(micro))
  for (i in seq_len(nrow(micro))){
    for (j in seq_len(ncol(micro))){
      val <- 2 * (micro[i,j] - centers[j]) / ranges[j]
      if (val == 1 || val == -1) expected[i,j] <- NA_real_ else expected[i,j] <- val
    }
  }

  expect_equal(U, expected)
})

test_that("get_latent_param returns expected values and errors", {
  res1 <- get_latent_param(LatentCase = "U_id_symmetric", LatentDist = "Unif", p = 2)
  expect_equal(res1$LatentParam[[1]], 1/12)
  expect_equal(res1$LatentCase, "U_id_symmetric")
  expect_equal(res1$LatentDist, "Unif")
  expect_equal(res1$TriangParam, 0)
  expect_equal(res1$BetaParam.a, 1)
  expect_equal(res1$BetaParam.b, 1)

  res2 <- get_latent_param(LatentCase = "U_id", LatentDist = "Triang", TriangParam = 0, p = 2)
  expect_equal(res2$LatentParam[[1]], 1/24)
  expect_equal(res2$LatentParam[[2]], 0)
  expect_equal(res2$LatentCase, "U_id")
  expect_equal(res2$LatentDist, "Triang")

  expect_error(get_latent_param(LatentCase = "General", LatentDist = "KDE", p = 2), "microdata must be provided")
  expect_error(get_latent_param(LatentCase = "U_id", LatentDist = "Unif", estimate.DistParam = TRUE, p = 2), "LatentCase to 'General'")
})

test_that("get_latent_param defaults LatentDist to Unif for symmetric case", {
  res <- get_latent_param(LatentCase = "U_id_symmetric", p = 3)
  expect_equal(res$LatentDist, "Unif")
  expect_equal(res$LatentParam[[1]], 1/12)
})

test_that("get_latent_param general case returns covariance matrix and mean matrix", {
  Umicro <- matrix(c(-0.5, 0.5, 0.1, -0.1, 0, 0), ncol = 2)
  res <- get_latent_param(LatentCase = "General", LatentDist = "Triang", TriangParam = c(0, 0.5), Umicro = Umicro)

  expect_equal(res$LatentCase, "General")
  expect_equal(res$LatentDist, "Triang")
  expect_length(res$LatentParam, 2)
  expect_true(is.matrix(res$LatentParam[[1]]))
  expect_equal(dim(res$LatentParam[[1]]), c(2, 2))
  expect_true(is.matrix(res$LatentParam[[2]]))
  expect_equal(dim(res$LatentParam[[2]]), c(2, 2))
  expect_equal(res$TriangParam, c(0, 0.5))
  expect_equal(res$BetaParam.a, 1)
  expect_equal(res$BetaParam.b, 1)
})

test_that("get_latent_param supports degenerated latent distribution", {
  expect_message(
    res <- get_latent_param(LatentCase = "U_id_symmetric", LatentDist = "Degenerated", p = 1),
    regexp="The latent variable, U, is a degenerated random variable"
  )
  
  expect_equal(res$LatentParam[[1]], 0)
  expect_equal(res$LatentDist, "Degenerated")
})

test_that("get_latent_param estimate.DistParam works for General case", {
  Umicro <- matrix(c(-0.4, 0.4, 0.2, -0.2), ncol = 2)
  res <- get_latent_param(LatentCase = "General", LatentDist = "KDE", Umicro = Umicro, estimate.DistParam = TRUE)

  expect_equal(res$LatentCase, "General")
  expect_equal(res$LatentDist, "KDE")
  expect_equal(res$TriangParam, c(NA, NA))
  expect_equal(res$BetaParam.a, c(NA, NA))
  expect_equal(res$BetaParam.b, c(NA, NA))
  expect_true(is.matrix(res$LatentParam[[1]]))
  expect_equal(dim(res$LatentParam[[1]]), c(2, 2))
})
 
test_that("get_latent_param requires microdata or p when needed", {
  expect_error(get_latent_param(LatentCase = "U_id", LatentDist = "Triang"), "The number of variables or the microdata must be provided")
  expect_error(get_latent_param(LatentCase = "General", LatentDist = "Triang", Umicro = NULL), "a distribution type must be provided|microdata must be provided")
})

test_that("get_latent_param rejects per-variable parameters outside General", {
  expect_error(get_latent_param(LatentCase = "U_id", LatentDist = "Triang", TriangParam = c(0, 0.5), p = 2), "LatentCase must be 'General'")
  expect_error(get_latent_param(LatentCase = "U_id", LatentDist = "Beta", BetaParam.a = c(1, 2), p = 2), "LatentCase must be 'General'")
})

test_that("get_latent_var handles different Seq formats correctly", {
  macro_LbUb <- data.frame(L1 = 1, U1 = 3, L2 = 0, U2 = 2)
  macro_AllLbAllUb <- data.frame(L1 = 1, L2 = 0, U1 = 3, U2 = 2)
  macro_CenRng <- data.frame(C1 = 2, R1 = 2, C2 = 1, R2 = 2)
  macro_AllCenAllRng <- data.frame(C1 = 2, C2 = 1, R1 = 2, R2 = 2)
  
  rownames(macro_LbUb) <- "G1"
  rownames(macro_AllLbAllUb) <- "G1"
  rownames(macro_CenRng) <- "G1"
  rownames(macro_AllCenAllRng) <- "G1"
  
  micro <- matrix(c(2, 1, 3, 2, 1, 0), ncol = 2, byrow = TRUE)
  agrby <- rep("G1", nrow(micro))

  U_LbUb <- get_latent_var(micro, macro_LbUb, agrby, "G1", Seq = "LbUb_VarbyVar")
  U_AllLbAllUb <- get_latent_var(micro, macro_AllLbAllUb, agrby, "G1", Seq = "AllLb_AllUb")
  U_CenRng <- get_latent_var(micro, macro_CenRng, agrby, "G1", Seq = "CenRng_VarbyVar")
  U_AllCenAllRng <- get_latent_var(micro, macro_AllCenAllRng, agrby, "G1", Seq = "AllCen_AllRng")

  # All formats should produce equivalent results
  expect_equal(U_LbUb, U_AllLbAllUb, tolerance = 1e-10)
  expect_equal(U_LbUb, U_CenRng, tolerance = 1e-10)
  expect_equal(U_LbUb, U_AllCenAllRng, tolerance = 1e-10)
})

test_that("get_latent_var handles multiple groups correctly", {
  macro <- data.frame(L1 = c(1, 5), U1 = c(3, 7), L2 = c(0, 2), U2 = c(2, 4))
  rownames(macro) <- c("G1", "G2")
  
  micro <- matrix(c(2, 1, 3, 2, 6, 3, 5, 2.5), ncol = 2, byrow = TRUE)
  agrby <- c("G1", "G1", "G2", "G2")

  U <- get_latent_var(micro, macro, agrby, rownames(macro), Seq = "LbUb_VarbyVar")

  expect_equal(nrow(U), nrow(micro))
  expect_equal(ncol(U), ncol(micro))
})

test_that("get_latent_var rejects invalid macrodata", {
  micro <- matrix(c(2, 1), ncol = 2)
  agrby <- "G1"
  
  expect_error(get_latent_var(micro, "invalid", agrby, "G1"), "macrodata must be a data.frame, matrix, or intData object")
})

test_that("get_latent_var handles edge values at boundaries", {
  # Test with micro values exactly at boundaries
  macro <- data.frame(L1 = 0, U1 = 2, L2 = 0, U2 = 2)
  rownames(macro) <- "G1"
  
  micro <- matrix(c(0, 0, 2, 2, 1, 1), ncol = 2, byrow = TRUE)
  agrby <- c("G1", "G1", "G1")

  U <- get_latent_var(micro, macro, agrby, "G1", Seq = "LbUb_VarbyVar")

  # Values at exact boundaries should be NA
  expect_true(is.na(U[1, 1]))  # U = -1
  expect_true(is.na(U[1, 2]))  # U = -1
  expect_true(is.na(U[2, 1]))  # U = 1
  expect_true(is.na(U[2, 2]))  # U = 1
  expect_false(is.na(U[3, 1])) # U = 0
})

test_that("get_latent_param with InvTri distribution", {
  res <- get_latent_param(LatentCase = "U_id_symmetric", LatentDist = "InvTri", p = 2)
  expect_equal(res$LatentParam[[1]], 1/8)
  expect_equal(res$LatentDist, "InvTri")
})

test_that("get_latent_param with TNorm distribution", {
  res <- get_latent_param(LatentCase = "U_id_symmetric", LatentDist = "TNorm", p = 1)
  expect_true(is.numeric(res$LatentParam[[1]]))
  expect_true(res$LatentParam[[1]] > 0 && res$LatentParam[[1]] < 0.1)
  expect_equal(res$LatentDist, "TNorm")
})

test_that("get_latent_param with U_id case and Beta distribution", {
  res_beta <- get_latent_param(LatentCase = "U_id", LatentDist = "Beta", BetaParam.a = 2, BetaParam.b = 3, p = 1)
  expect_true(is.numeric(res_beta$LatentParam[[1]]))
  expect_true(is.numeric(res_beta$LatentParam[[2]]))
})

test_that("get_latent_param with multidimensional General case", {
  Umicro <- matrix(c(-0.5, 0.5, 0.2, 0.1, -0.2, 0, 0.3, 0.1, -0.1), ncol = 3)
  res <- get_latent_param(LatentCase = "General", LatentDist = "KDE", Umicro = Umicro)
  
  expect_equal(res$LatentCase, "General")
  expect_equal(dim(res$LatentParam[[1]]), c(3, 3))  # E_UU matrix
  expect_equal(dim(res$LatentParam[[2]]), c(3, 3))  # Psi matrix
})

test_that("get_latent_param General case with mixed distributions", {
  Umicro <- matrix(c(-0.5, 0.5, 0.2, 0.1, -0.2, 0), ncol = 2)
  res <- get_latent_param(LatentCase = "General", 
                         LatentDist = c("Beta", "KDE"),
                         BetaParam.a = c(1.5, NA),
                         BetaParam.b = c(2, NA),
                         Umicro = Umicro)
  
  expect_equal(res$LatentDist, c("Beta", "KDE"))
  expect_equal(res$BetaParam.a[1], 1.5)
})

test_that("get_latent_param with Triang distribution different modes", {
  res0 <- get_latent_param(LatentCase = "U_id", LatentDist = "Triang", TriangParam = 0, p = 1)
  res_pos <- get_latent_param(LatentCase = "U_id", LatentDist = "Triang", TriangParam = 0.5, p = 1)
  res_neg <- get_latent_param(LatentCase = "U_id", LatentDist = "Triang", TriangParam = -0.5, p = 1)
  
  expect_equal(res0$LatentParam[[2]], 0)  # mean should be 0 for mode=0
  expect_equal(res_pos$LatentParam[[2]], 0.5/3)  # mean = mode/3
  expect_equal(res_neg$LatentParam[[2]], -0.5/3)
})

test_that("get_latent_param General case with per-variable Beta parameters", {
  Umicro <- matrix(c(-0.4, 0.4, 0.2, -0.2, 0.1, 0.3), ncol = 3)
  res <- get_latent_param(LatentCase = "General",
                         LatentDist = c("Beta", "Beta", "KDE"),
                         BetaParam.a = c(1.5, 2, NA),
                         BetaParam.b = c(2, 3, NA),
                         Umicro = Umicro)
  
  expect_equal(res$BetaParam.a[1], 1.5)
  expect_equal(res$BetaParam.a[2], 2)
  expect_equal(dim(res$LatentParam[[1]]), c(3, 3))
})

test_that("get_latent_param with per-variable Triang parameters", {
  Umicro <- matrix(c(-0.5, 0.5, 0.1, -0.1), ncol = 2)
  res <- get_latent_param(LatentCase = "General",
                         LatentDist = "Triang",
                         TriangParam = c(0, 0.3),
                         Umicro = Umicro)
  
  expect_equal(res$TriangParam, c(0, 0.3))
  expect_equal(dim(res$LatentParam[[1]]), c(2, 2))
})

test_that("meanU returns scalar for single distribution", {
  m_unif <- AIDA:::meanU(LatentDist = "Unif", p = 1)
  expect_equal(m_unif, 0)
  
  m_triang <- AIDA:::meanU(LatentDist = "Triang", TriangParam = 0.5, p = 1)
  expect_equal(m_triang, 0.5/3)
})

test_that("meanU returns matrix for multiple distributions", {
  m_multi <- AIDA:::meanU(LatentDist = c("Unif", "Triang", "Beta"),
                   TriangParam = c(NA, 0.2, NA),
                   BetaParam.a = c(NA, NA, 2),
                   BetaParam.b = c(NA, NA, 3),
                   p = 3)
  
  expect_true(is.matrix(m_multi))
  expect_equal(dim(m_multi), c(3, 3))
  expect_equal(diag(m_multi)[1], 0)  # Unif mean
  expect_equal(diag(m_multi)[2], 0.2/3)  # Triang mean
})

test_that("meanU with KDE distribution", {
  Umicro <- matrix(c(-0.5, 0.5, 0.1, -0.1, 0, 0.2), ncol = 2)
  m_kde <- AIDA:::meanU(LatentDist = c("KDE", "Unif"), Umicro = Umicro, p = 2)
  
  expect_true(is.matrix(m_kde))
  expect_equal(m_kde[1, 1], mean(Umicro[, 1], na.rm = TRUE))
  expect_equal(m_kde[2, 2], 0)
})

test_that("meanU2 returns scalar for single distribution", {
  m2_unif <- AIDA:::meanU2(LatentDist = "Unif", p = 1)
  expect_equal(m2_unif, 1/3)
  
  m2_triang <- AIDA:::meanU2(LatentDist = "Triang", TriangParam = 0, p = 1)
  expect_equal(m2_triang, 1/6)
})

test_that("meanU2 returns matrix for multiple distributions", {
  m2_multi <- AIDA:::meanU2(LatentDist = c("Unif", "Triang", "Beta"),
                     TriangParam = c(NA, 0, NA),
                     BetaParam.a = c(NA, NA, 2),
                     BetaParam.b = c(NA, NA, 3),
                     p = 3)
  
  expect_true(is.matrix(m2_multi))
  expect_equal(dim(m2_multi), c(3, 3))
  expect_equal(diag(m2_multi)[1], 1/3)  # Unif
  expect_equal(diag(m2_multi)[2], 1/6)  # Triang with mode=0
})

test_that("meanU2 with KDE distribution", {
  Umicro <- matrix(c(-0.5, 0.5, 0.1, -0.1, 0, 0.2), ncol = 2)
  m2_kde <- AIDA:::meanU2(LatentDist = c("KDE", "Unif"), Umicro = Umicro, p = 2)
  
  expect_true(is.matrix(m2_kde))
  expect_equal(m2_kde[1, 1], mean(Umicro[, 1]^2, na.rm = TRUE))
  expect_equal(m2_kde[2, 2], 1/3)
})

test_that("cal.E.UU with Beta-Beta distributions", {
  e_uu <- AIDA:::cal.E.UU(LatentDist = "Beta",
                   BetaParam.a = c(1, 2),
                   BetaParam.b = c(1, 2),
                   p = 2)
  
  expect_true(is.matrix(e_uu))
  expect_equal(dim(e_uu), c(2, 2))
  expect_equal(e_uu[1, 2], e_uu[2, 1])  # Symmetry
})

test_that("cal.E.UU with Triang-Triang distributions", {
  e_uu <- AIDA:::cal.E.UU(LatentDist = "Triang",
                   TriangParam = c(0, 0.2),
                   p = 2)
  
  expect_true(is.matrix(e_uu))
  expect_equal(dim(e_uu), c(2, 2))
  expect_equal(e_uu[1, 2], e_uu[2, 1])
})

test_that("cal.E.UU with KDE distributions", {
  Umicro <- matrix(c(-0.5, 0.5, 0.1, -0.1, 0, 0.2), ncol = 2)
  e_uu <- AIDA:::cal.E.UU(LatentDist = "KDE", Umicro = Umicro, p = 2)
  
  expect_true(is.matrix(e_uu))
  expect_equal(dim(e_uu), c(2, 2))
  expect_equal(e_uu[1, 2], e_uu[2, 1])  # Symmetry
  expect_true(all(is.finite(e_uu)))
})

test_that("CalE.beta.beta computes correctly", {
  calE <- AIDA:::CalE.beta.beta(a1 = 1, b1 = 1, a2 = 1, b2 = 1)
  
  expect_true(is.numeric(calE))
  expect_true(is.finite(calE))
  expect_true(calE >= 0)
})

test_that("CalE.beta.kde computes correctly", {
  micro <- c(-0.5, -0.3, -0.1, 0.1, 0.3, 0.5)
  calE <- AIDA:::CalE.beta.kde(micro, a1 = 1, b1 = 1)
  
  expect_true(is.numeric(calE))
  expect_true(is.finite(calE))
})

test_that("CalE.kde.kde computes correctly", {
  micro1 <- c(-0.5, -0.3, -0.1, 0.1, 0.3, 0.5)
  micro2 <- c(-0.4, -0.2, 0, 0.2, 0.4, 0.6)
  calE <- AIDA:::CalE.kde.kde(micro1, micro2)
  
  expect_true(is.numeric(calE))
  expect_true(is.finite(calE))
  expect_true(calE >= 0)
})

test_that("CalE.triang.triang computes correctly", {
  calE_00 <- AIDA:::CalE.triang.triang(mo1 = 0, mo2 = 0)
  expect_equal(calE_00, 1/6)
  
  calE_diff <- AIDA:::CalE.triang.triang(mo1 = 0.3, mo2 = -0.2)
  expect_true(is.numeric(calE_diff))
  expect_true(is.finite(calE_diff))
})

test_that("get_latent_param with estimate.DistParam for Triang", {
  Umicro <- matrix(c(-0.3, 0.3, 0.2, 0.3, -0.2, 0.4), ncol = 2)
  res <- get_latent_param(LatentCase = "General",
                         LatentDist = c("Triang", "Triang"),
                         Umicro = Umicro,
                         estimate.DistParam = TRUE)
  
  expect_equal(res$LatentCase, "General")
  expect_true(is.numeric(res$TriangParam))
  expect_equal(length(res$TriangParam), 2)
})

test_that("get_latent_param with estimate.DistParam for Beta", {
  Umicro <- matrix(c(-0.3, 0.3, 0.2, 0.3, -0.2, 0.4), ncol = 2)
  res <- get_latent_param(LatentCase = "General",
                         LatentDist = c("Beta","Beta"),
                         Umicro = Umicro,
                         estimate.DistParam = TRUE)
  
  expect_equal(res$LatentCase, "General")
  expect_true(all(is.numeric(res$BetaParam.a)))
  expect_true(all(is.numeric(res$BetaParam.b)))
  expect_equal(length(res$BetaParam.a), 2)
  expect_equal(length(res$BetaParam.b), 2)
})

test_that("get_latent_param error when estimate.DistParam=TRUE without Umicro", {
  expect_error(get_latent_param(LatentCase = "General",
                               LatentDist = "Beta",
                               estimate.DistParam = TRUE,
                               p = 2),
              "microdata must be provided")
})

test_that("get_latent_param error for invalid Seq format in get_latent_var", {
  macro <- data.frame(L1 = 1, U1 = 3, L2 = 0, U2 = 2)
  rownames(macro) <- "G1"
  micro <- matrix(c(2, 1), ncol = 2)
  
  expect_error(get_latent_var(micro, macro, "G1", "G1", Seq = "Invalid"),
              "should be one of")
})

test_that("get_latent_var preserves data types and dimensions", {
  macro <- data.frame(L1 = 1, U1 = 3, L2 = 0, U2 = 2)
  rownames(macro) <- "G1"
  micro <- matrix(c(2, 1, 3, 2, 1, 0.5), ncol = 2, byrow = TRUE)
  agrby <- rep("G1", nrow(micro))

  U <- get_latent_var(micro, macro, agrby, "G1", Seq = "LbUb_VarbyVar")

  expect_equal(class(U), c("matrix","array"))
  expect_equal(nrow(U), nrow(micro))
  expect_equal(ncol(U), ncol(micro))
  expect_true(is.numeric(U))
})

test_that("get_latent_var with single variable", {
  macro <- data.frame(L = 0, U = 10)
  rownames(macro) <- "G1"
  micro <- matrix(c(2, 5, 8), ncol = 1)
  agrby <- c("G1", "G1", "G1")

  U <- get_latent_var(micro, macro, agrby, "G1", Seq = "LbUb_VarbyVar")

  expect_equal(nrow(U), 3)
  expect_equal(ncol(U), 1)
  expected <- 2 * (micro - 5) / 10
  expect_equal(U, expected)
})

test_that("get_latent_var with many variables", {
  macro <- data.frame(L1 = 1, U1 = 3, L2 = 0, U2 = 2, L3 = -1, U3 = 1, L4 = 10, U4 = 20)
  rownames(macro) <- "G1"
  micro <- matrix(runif(20, 0, 20), ncol = 4)
  agrby <- rep("G1", nrow(micro))

  U <- get_latent_var(micro, macro, agrby, "G1", Seq = "LbUb_VarbyVar")

  expect_equal(nrow(U), nrow(micro))
  expect_equal(ncol(U), 4)
})
