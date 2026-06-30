library(testthat)

with_null_device <- function(code) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  force(code)
}

make_shapley_matrix <- function(n = 3, p = 3) {
  mat <- matrix(seq(1, n * p), nrow = n, ncol = p, byrow = TRUE)
  rownames(mat) <- paste0("obs", seq_len(n))
  colnames(mat) <- paste0("V", seq_len(p))
  mat
}

test_that("plot_int_Shapley_inter draws without error", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("cowplot")
  x <- matrix(c(1,2,3, -1,0,2, 0.5,0.5,-0.2), nrow = 3, byrow = TRUE)
  rownames(x) <- colnames(x) <- c("A","B","C")
  expect_s3_class(plot_int_Shapley_inter(x), "ggplot")
  expect_error(with_null_device(plot_int_Shapley_inter(x)), NA)
})

test_that("plot_bar_int_Shapley returns ggplot and errors on bad cutoff", {
  testthat::skip_if_not_installed("ggplot2")
  shp <- make_shapley_matrix(4, 3)
  res <- plot_bar_int_Shapley(shp, palette = rainbow(ncol(shp)))
  expect_s3_class(res, "ggplot")

  expect_error(plot_bar_int_Shapley(shp, cutoff_value = "bad", palette = rainbow(ncol(shp))), "cutoff_value must be numeric")
})

test_that("plot_tile_int_Shapley and plot_radar_int_Shapley run", {
  testthat::skip_if_not_installed("ggplot2")
  shp <- make_shapley_matrix(4, 3)
  expect_s3_class(plot_tile_int_Shapley(shp), "ggplot")

  # radarplot uses base graphics and fmsb
  testthat::skip_if_not_installed("fmsb")
  expect_error(with_null_device(plot_radar_int_Shapley(shp)), NA)
})

test_that("plot_beeswarm_int_Shapley and barplot_decomp produce plots when deps present", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggbeeswarm")
  testthat::skip_if_not_installed("plotly")

  shp <- make_shapley_matrix(5, 3)
  color_class <- rep(c("Regular","Outlier"), length.out = nrow(shp))
  res <- plot_beeswarm_int_Shapley(shp, color_class = color_class, ggplotly = FALSE)
  expect_s3_class(res, "ggplot")

  # plot_bar_int_Shapley_decomp expects a list of matrices
  decomp <- list(obs1 = matrix(1:6, nrow = 2, byrow = TRUE), obs2 = matrix(2:7, nrow = 2, byrow = TRUE))
  colnames(decomp[[1]]) <- colnames(decomp[[2]]) <- c("Centers","Ranges","CentersRanges")[seq_len(ncol(decomp[[1]]))]
  rownames(decomp[[1]]) <- rownames(decomp[[2]]) <- c("X","Y")
  res2 <- plot_bar_int_Shapley_decomp(decomp, plot_IMah = FALSE, palette = rainbow(3))
  expect_s3_class(res2, "ggplot")
})

test_that("plot_beeswarm_int_Shapley ggplotly returns plotly object", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("plotly")
  shp <- make_shapley_matrix(5, 3)
  color_class <- rep(c("Regular","Outlier"), length.out = nrow(shp))
  res_plotly <- plot_beeswarm_int_Shapley(shp, color_class = color_class, ggplotly = TRUE)
  expect_s3_class(res_plotly, "plotly")
})

test_that("plots error on clearly invalid inputs", {
  testthat::skip_if_not_installed("ggplot2")
  expect_error(plot_int_Shapley_inter(1:5))
  expect_error(plot_bar_int_Shapley(list(1,2)))
  expect_error(plot_tile_int_Shapley(1:5))
  expect_error(plot_radar_int_Shapley(list(1,2)))
  expect_error(plot_beeswarm_int_Shapley(1:5))
  expect_error(plot_bar_int_Shapley_decomp(matrix(1:9, 3)))
})

test_that("plot_bar_int_Shapley handles cutoff values and plot_IMah FALSE", {
  testthat::skip_if_not_installed("ggplot2")
  shp <- make_shapley_matrix(4, 3)
  p <- plot_bar_int_Shapley(shp, plot_IMah = FALSE, cutoff_value = c(1, 2), cutoff_label = c("L1", "L2"), rotate_x = FALSE, palette = rainbow(ncol(shp)))
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))))
  expect_error(plot_bar_int_Shapley(shp, cutoff_value = c(1, 2), cutoff_label = "bad", palette = rainbow(ncol(shp))), "same length")
})

test_that("plot_tile_int_Shapley highlights outliers and show_values", {
  testthat::skip_if_not_installed("ggplot2")
  shp <- make_shapley_matrix(4, 3)
  outliers <- list(outliers_names = c("obs1", "obs3"))
  p <- plot_tile_int_Shapley(shp, outliers = outliers, show_values = TRUE, sort.var = TRUE, sort.obs = TRUE)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)
})

test_that("plot_radar_int_Shapley draws without error for sorted observations", {
  testthat::skip_if_not_installed("fmsb")
  shp <- make_shapley_matrix(4, 3)
  expect_error(with_null_device(plot_radar_int_Shapley(shp, palette = c("red", "blue", "green", "purple"), sort.obs = TRUE)), NA)
})

test_that("plot_beeswarm_int_Shapley supports ggplotly = FALSE and shape labels", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggbeeswarm")
  testthat::skip_if_not_installed("ggrepel")
  shp <- make_shapley_matrix(5, 3)
  color_class <- rep(c("Regular", "Outlier"), length.out = nrow(shp))
  shape_class <- rep(c("A", "B"), length.out = nrow(shp))
  p <- plot_beeswarm_int_Shapley(shp, color_class = color_class, shape_class = shape_class,
                            color_label = "Class", shape_label = "Shape",
                            ggplotly = FALSE, label_obs = "obs1")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$colour, "Class")
  expect_equal(p$labels$shape, "Shape")
})

test_that("plot_bar_int_Shapley_decomp can plot IMah cutoff lines", {
  testthat::skip_if_not_installed("ggplot2")
  decomp <- list(obs1 = matrix(1:6, nrow = 2, byrow = TRUE), obs2 = matrix(2:7, nrow = 2, byrow = TRUE))
  colnames(decomp[[1]]) <- colnames(decomp[[2]]) <- c("Centers", "Ranges", "CentersRanges")[seq_len(ncol(decomp[[1]]))]
  rownames(decomp[[1]]) <- rownames(decomp[[2]]) <- c("X", "Y")
  p <- plot_bar_int_Shapley_decomp(decomp, plot_IMah = TRUE, palette = rainbow(3))
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))))
})

test_that("ggplot outputs contain layers and labels", {
  testthat::skip_if_not_installed("ggplot2")
  shp <- make_shapley_matrix(4, 3)
  p <- plot_bar_int_Shapley(shp, palette = rainbow(ncol(shp)))
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1)

  t <- plot_tile_int_Shapley(shp)
  expect_s3_class(t, "ggplot")
  expect_true(length(t$layers) >= 1)
})
