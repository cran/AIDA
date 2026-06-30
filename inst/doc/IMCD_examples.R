## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----results='hide', message = FALSE, warning=FALSE---------------------------
library(AIDA)

## -----------------------------------------------------------------------------
data(intCars)
cars_microdata <- intCars$microdata
cars_int <- intCars$intData

## -----------------------------------------------------------------------------
cars_IMCD <- IMCD(cars_int, m = floor(0.75*cars_int@NObs), cutoff = "farness", cutoff_lvl = 0.9)
cars_outliers <- int_outliers(cars_IMCD$robust_dist, cutoff = "farness", cutoff_lvl = 0.9)
cars_outliers$outliers_names

## -----------------------------------------------------------------------------
cars_outliers_colors <- rep('gray50', cars_int@NObs)
names(cars_outliers_colors) <- rownames(cars_int)
cars_outliers_colors[cars_outliers$outliers_names] <- 'red'

plot_pairs_int(cars_int, palette = cars_outliers_colors, type = "rectangles", 
                    corr = cov2cor(cars_IMCD$cov_IMCD), labels = colnames(cars_int),
                    is_outlier = cars_outliers$is_outlier, gap = 0)

## ----eval = requireNamespace("ggrepel", quietly = TRUE) && requireNamespace("robustbase", quietly = TRUE)----
# Classical distances and outliers
cars_class_dist <- IMah_dist(cars_int, z = rep(1,cars_int@NObs))
cars_class_outliers <- int_outliers(cars_class_dist, cutoff = "adjbox", cutoff_lvl = 1.5)

cars_is_outliers <- as.character(cars_outliers$is_outlier)
cars_is_outliers[cars_outliers$is_outlier] <- "Outlier"
cars_is_outliers[!cars_outliers$is_outlier] <- "Inlier"

plot_dist_dist(cars_class_dist, cars_class_outliers$cutoff_value[[2]], class_cutoff_label = "1.5 Adjbox",
                cars_IMCD$robust_dist, cars_outliers$cutoff_value, rob_cutoff_label = "0.9 Farness",
                color_class = cars_is_outliers, ggplotly = FALSE, shape_class = cars_microdata$class, 
                shape_label = "Class", palette = c("gray50","red"), 
                label_obs = c(cars_outliers$outliers_names, "Bmwserie7"))

## -----------------------------------------------------------------------------
data(spotify_tracks)
spotify_int <- spotify_tracks$intData_trimmed

## -----------------------------------------------------------------------------
spotify_IMCD <- IMCD(spotify_int, m = round(0.75*nrow(spotify_int)), 
                      cutoff = "farness", cutoff_lvl = 0.95)

# Strong outliers
spotify_outliers <- int_outliers(spotify_IMCD$robust_dist, cutoff = "farness", cutoff_lvl = 0.95)
spotify_outliers$outliers_names

# Mild outliers
spotify_outliers_2 <- int_outliers(spotify_IMCD$robust_dist, cutoff = "farness", cutoff_lvl = 0.9)
spotify_outliers_2$outliers_names[!spotify_outliers_2$outliers_names%in%spotify_outliers$outliers_names]

## ----eval = requireNamespace("corrplot", quietly = TRUE), fig.width=9---------
# Compute correlation matrix from the robust covariance matrix
spotify_corr <- cov2cor(spotify_IMCD$cov_IMCD)

colfunc <- colorRampPalette(c("deepskyblue4", "white", "red4"))
corrplot::corrplot(
  spotify_corr,
  method = "color",
  type = "upper",
  diag = FALSE,
  col = colfunc(200),
  tl.col = "black",
  tl.srt = 45,
  tl.offset = 1,
  tl.cex = 1.2,
  outline = TRUE,
  addCoef.col = "black",
  number.cex = 1
)

## ----eval = requireNamespace("ggrepel", quietly = TRUE) && requireNamespace("robustbase", quietly = TRUE)----
# Classical distances and outliers
spotify_class_dist <- IMah_dist(spotify_int, z = rep(1,spotify_int@NObs))
spotify_class_outliers <- int_outliers(spotify_class_dist, cutoff = "adjbox")

spotify_is_outliers <- as.character(spotify_outliers$is_outlier)
spotify_is_outliers[!spotify_outliers_2$is_outlier] <- "Regular"
spotify_is_outliers[spotify_outliers_2$is_outlier] <- "Mild Outlier"
spotify_is_outliers[spotify_outliers$is_outlier] <- "Extreme Outlier"
palette_outliers <- c(
  "Regular"         = "gray50",
  "Mild Outlier"    = "darkorange",
  "Extreme Outlier" = "red"
)

plot_dist_dist(spotify_class_dist, spotify_class_outliers$cutoff_value[[2]], 
                "1.5 Adjbox", spotify_IMCD$robust_dist, 
                c(spotify_outliers_2$cutoff_value, spotify_outliers$cutoff_value), 
                c("0.90 Farness", "0.95 Farness"), color_class = spotify_is_outliers, 
                color_label = "Outlier Status", palette = palette_outliers, ggplotly = FALSE, 
                label_obs = spotify_outliers_2$outliers_names)

