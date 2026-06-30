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
cars_IMCD <- IMCD(cars_int, m=floor(0.75*cars_int@NObs), cutoff = "farness", cutoff_lvl = 0.9)
cars_outliers <- int_outliers(cars_IMCD$robust_dist, p = cars_int@NIVar,
                                    cutoff = "farness", cutoff_lvl = 0.9)
cars_outliers$outliers_names

## ----fig.width=7, fig.height=4------------------------------------------------
cars_is_outliers <- as.character(cars_outliers$is_outlier)
cars_is_outliers[cars_outliers$is_outlier] <- "Outlier"
cars_is_outliers[!cars_outliers$is_outlier] <- "Inlier"

plot_interval_dist(
  dist = cars_IMCD$robust_dist,
  cutoff = cars_outliers$cutoff_value,
  cutoff_label = "0.9 Farness",
  obs_names = rownames(cars_int),
  color_class = cars_is_outliers,
  palette =c("gray50","dodgerblue"),
  shape_class = cars_microdata$class,
  shape_label = "Class",
  sort.obs = TRUE
)

## -----------------------------------------------------------------------------
cars_shapley <- int_Shapley(cars_int, mean_c = cars_IMCD$mean_IMCD_c, 
                            mean_r = cars_IMCD$mean_IMCD_r, cov = cars_IMCD$cov_IMCD)

cars_shapley2 <- int_Shapley(cars_int)

## ----eval = requireNamespace("scales", quietly = TRUE)------------------------
plot_bar_int_Shapley(cars_shapley[c(cars_outliers$outliers_names,"Bmwserie7"),], 
                    cutoff_value = cars_outliers$cutoff_value, 
                    cutoff_label = "0.9 Farness Cutoff", 
                    palette = scales::hue_pal()(4))

## ----eval = requireNamespace("ggrepel", quietly = TRUE)-----------------------
plot_beeswarm_int_Shapley(cars_shapley, cars_is_outliers, color_label = NULL, 
                      shape_class = cars_microdata$class, shape_label = "Class",
                      palette = c("gray50","dodgerblue"), ggplotly = FALSE, 
                      label_obs = c(cars_outliers$outliers_names), rotate_x = FALSE)

## ----fig.width=9.5, fig.height=4----------------------------------------------
plot_tile_int_Shapley(cars_shapley, abbrev.var = 15, sort.obs = TRUE)

## ----fig.width=9.5------------------------------------------------------------
outliers_colors <- rep('gray50', cars_int@NObs)
names(outliers_colors) <- rownames(cars_int)
outliers_colors[cars_outliers$outliers_names] = 'dodgerblue'

plot_radar_int_Shapley(cars_shapley,outliers_colors)

## ----fig.width=7--------------------------------------------------------------
cars_shapley_inter <- int_Shapley_interaction(cars_int, mean_c = cars_IMCD$mean_IMCD_c,
                                              mean_r = cars_IMCD$mean_IMCD_r, 
                                              cov = cars_IMCD$cov_IMCD)

plot_int_Shapley_inter(cars_shapley_inter[["Ferrari"]], abbrev = 15, title = "Ferrari")

## ----eval = requireNamespace("scales", quietly = TRUE), fig.width=9.5, fig.height=5----
cars_shapley_decomp <- int_Shapley_decomp(cars_int, mean_c = cars_IMCD$mean_IMCD_c,
                                          mean_r = cars_IMCD$mean_IMCD_r, cov = cars_IMCD$cov_IMCD)

plot_bar_int_Shapley_decomp(cars_shapley_decomp[c(cars_outliers$outliers_names,"Bmwserie7")],
                          rotate_x = FALSE, palette = scales::hue_pal()(4))

## -----------------------------------------------------------------------------
data(spotify_tracks)
spotify_int <- spotify_tracks$intData_trimmed

## -----------------------------------------------------------------------------
spotify_IMCD <- IMCD(spotify_int, m=round(0.75*nrow(spotify_int)), 
                              cutoff="farness", cutoff_lvl = 0.95)

# Strong outliers
spotify_outliers <- int_outliers(spotify_IMCD$robust_dist,cutoff="farness", cutoff_lvl = 0.95)
spotify_outliers$outliers_names

# Mild outliers
spotify_outliers_2 <- int_outliers(spotify_IMCD$robust_dist,cutoff="farness", cutoff_lvl = 0.9)
spotify_outliers_2$outliers_names[!spotify_outliers_2$outliers_names%in%spotify_outliers$outliers_names]

## ----eval = requireNamespace("ggrepel", quietly = TRUE), fig.width=7, fig.height=4----
spotify_is_outliers <- as.character(spotify_outliers$is_outlier)
spotify_is_outliers[!spotify_outliers_2$is_outlier] <- "Regular"
spotify_is_outliers[spotify_outliers_2$is_outlier] <- "Mild Outlier"
spotify_is_outliers[spotify_outliers$is_outlier] <- "Extreme Outlier"
palette_outliers <- c(
  "Regular"        = "gray50",
  "Mild Outlier"   = 'forestgreen', #"darkorange",
  "Extreme Outlier"= 'dodgerblue' #"red"
)

plot_interval_dist(
  spotify_IMCD$robust_dist,
  c(spotify_outliers$cutoff_value,spotify_outliers_2$cutoff_value),
  c("0.95 Farness", "0.90 Farness"),
  sort.obs = FALSE,
  color_class = spotify_is_outliers, 
  color_label = "Outlier Status", 
  palette = palette_outliers,
  label_obs = spotify_outliers_2$outliers_names)

## ----eval = requireNamespace("ggrepel", quietly = TRUE) && requireNamespace("RColorBrewer", quietly = TRUE)----
spotify_Shapley <- int_Shapley(spotify_int, mean_c = spotify_IMCD$mean_IMCD_c, 
                              mean_r = spotify_IMCD$mean_IMCD_r, cov = spotify_IMCD$cov_IMCD)

high_dist_12 <- names(spotify_IMCD$robust_dist[order(spotify_IMCD$robust_dist, decreasing = TRUE)[1:12]])

plot_bar_int_Shapley(spotify_Shapley[high_dist_12,], 
                    cutoff_value = c(spotify_outliers$cutoff_value, spotify_outliers_2$cutoff_value),
                    cutoff_label = c("0.95 Farness Cutoff", "0.90 Farness Cutoff"),
                    sort.obs = TRUE, abbrev.obs = 20)

## -----------------------------------------------------------------------------
plot_beeswarm_int_Shapley(spotify_Shapley, spotify_is_outliers, "Outlier Status", 
                      palette = palette_outliers, ggplotly = FALSE, 
                      label_obs = c("sleep","classical"))

## ----fig.width=7--------------------------------------------------------------
plot_tile_int_Shapley(spotify_Shapley[high_dist_12,], sort.obs = TRUE, abbrev.var = 15)

## ----fig.width=7--------------------------------------------------------------
spotify_shapley_inter <- int_Shapley_interaction(spotify_int, mean_c = spotify_IMCD$mean_IMCD_c,
                                                  mean_r = spotify_IMCD$mean_IMCD_r, 
                                                  cov = spotify_IMCD$cov_IMCD)

plot_int_Shapley_inter(spotify_shapley_inter[["grindcore"]], abbrev = 20, title = "grindcore")

## ----eval = requireNamespace("RColorBrewer", quietly = TRUE), fig.width=9.5, fig.height=5----
spotify_shapley_decomp <- int_Shapley_decomp(spotify_int, mean_c = spotify_IMCD$mean_IMCD_c, 
                                              mean_r = spotify_IMCD$mean_IMCD_r, 
                                              cov = spotify_IMCD$cov_IMCD)

plot_bar_int_Shapley_decomp(spotify_shapley_decomp[spotify_outliers_2$outliers_names])

