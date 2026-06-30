#' Plot Shapley interaction indices
#'
#' @param x A \eqn{p \times p} matrix containing the Shapley interaction indices of a single observation.
#' @param abbrev Integer. If \code{abbrev.var} \eqn{> 0}, variable names are abbreviated using abbreviate with \code{minlenght = abrev}.
#' @param title Character. Title of the plot.
#' @param legend Logical. If TRUE (default), a legend is plotted.
#' @param text_size Integer. Size of the text in the plot 
#' @return Returns a figure consisting of two panels. The right panel shows the Shapley values, and the left panel the Shapley interaction indices.
#' @references Adapted from package \code{ShapleyOutlier} (\url{https://CRAN.R-project.org/package=ShapleyOutlier}).
#' @export
#'
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Estimate the mean and covariance matrix
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
#' 
#' # Compute Shapley interaction indices
#' credit_card_shap_inter <- int_Shapley_interaction(credit_card_int, 
#'                                                   mean_c = credit_card_IMCD$mean_IMCD_c, 
#'                                                   mean_r = credit_card_IMCD$mean_IMCD_r, 
#'                                                   cov = credit_card_IMCD$cov_IMCD)
#'
#' # Plot Shapley interaction for 1st observation
#' plot_int_Shapley_inter(credit_card_shap_inter[[1]])
plot_int_Shapley_inter <- function(x, 
                                  abbrev = 10, 
                                  title = NULL, 
                                  legend = TRUE, 
                                  text_size = 22){

  if(!is.matrix(x)) stop("`x` must be a matrix.")

  rowname <- name <- value <- Observation <- Value <- NULL

  # Abbreviate names
  if(abbrev > 0 && !is.null(rownames(x)) && !is.null(colnames(x))){
    rownames(x) <- abbreviate(rownames(x), minlength = abbrev)
    colnames(x) <- abbreviate(colnames(x), minlength = abbrev)
  }

  # Default names if missing
  if(is.null(rownames(x))){
    rownames(x) <- paste0("X", formatC(seq_len(nrow(x)), width = nchar(nrow(x)), flag = "0"))
  }
  if(is.null(colnames(x))){
    colnames(x) <- paste0("X", formatC(seq_len(ncol(x)), width = nchar(ncol(x)), flag = "0"))
  }

  # Reshape
  plot_PHI <- as.data.frame(as.table(x))
  colnames(plot_PHI) <- c("rowname", "name", "value")

  plot_PHI$value <- as.numeric(plot_PHI$value)
  plot_PHI$name <- factor(plot_PHI$name, levels = colnames(x))
  plot_PHI$rowname <- factor(plot_PHI$rowname, levels = rownames(x))

  # Heatmap
  pp <- ggplot(plot_PHI, aes(x = name, y = rowname, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "deepskyblue4", mid = "white",
                         high = "red4", midpoint = 0) + 
    theme_light() +
    theme(text = element_text(size = text_size),
                    axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1),
                    legend.text = element_text(size = text_size * 0.8),
                    legend.title = element_text(size = text_size * 0.8),
                    legend.position = "top",
                    legend.key.width = grid::unit(2, "cm"),
                    title = element_text(size = text_size * 0.8)) +
    labs(x = NULL, y = NULL, fill = "Shapley Inter.", title = title) + 
    scale_y_discrete(limits = rev(rownames(x)), expand = c(0, 0))

  if(!legend){
    pp <- pp + theme(legend.position = "none")
  }

  # Side barplot
  shval_df <- data.frame(
    Observation = factor(rownames(x), levels = rev(rownames(x))),
    Value = rowSums(x)
  )

  mp <- ggplot(shval_df, aes(x = Value, y = Observation)) +
    geom_col(width = 0.95) +
    theme_light() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = text_size),
          axis.title.x = element_text(size = text_size * 0.8)) +
    labs(x = "Shapley V.") +
    scale_y_discrete(limits = rev(rownames(x)), expand = c(0, 0))

  # Combine
  aligned_plots <- cowplot::align_plots(mp, pp, align = "h", axis = "tb")
  cowplot::plot_grid(aligned_plots[[2]], aligned_plots[[1]], ncol = 2, rel_widths = c(3, 1))
}

#' Barplot of Shapley values for Interval-valued Data
#'
#' @param x A \eqn{n \times p} matrix containing the Shapley values of \eqn{n} observations and \eqn{p} variables.
#' @param cutoff_value Numeric. The cutoff value used for detecting outliers. If \code{cutoff_value} is not \code{NULL} (default), the cutoff value is included in the plot.
#' @param cutoff_label Character. Label for the cutoff value line in the plot.
#' @param palette A vector with colors for each variable. If \code{palette} is \code{NULL} (default), the colors are generated using \code{RColorBrewer}.
#' @param abbrev.var Integer. If \code{abbrev.var} \eqn{> 0}, column names are abbreviated using abbreviate with \code{minlenght = abrev.var}.
#' @param abbrev.obs Integer. If \code{abbrev.obs} \eqn{> 0}, row names are abbreviated using abbreviate with \code{minlenght = abrev.obs}.
#' @param sort.obs Logical. If \code{TRUE} (default), observations are sorted according to their squared (robust) Interval-Mahalanobis distance.
#' @param plot_IMah Logical. If \code{TRUE} (default), the squared (robust) Interval-Mahalanobis distance will be included in the plot.
#' @param IMah_label Character. Label for the Interval-Mahalanobis distance in the plot legend. Default is "Robust \eqn{d_\mathrm{IMah}^2(\boldsymbol{x})}".
#' @param rotate_x Logical. If \code{TRUE} (default), the x-axis labels are rotated.
#'
#' @return Returns a barplot that displays the Shapley values (\code{\link{int_Shapley}}) for each observation and optionally (\code{plot_IMah = TRUE})
#' includes the squared (robust) Interval-Mahalanobis distance (\code{\link{IMah_dist}}) (black bar) and the corresponding outlier detection cut-off value (dotted line).
#' @references Adapted from package \code{ShapleyOutlier} (\url{https://CRAN.R-project.org/package=ShapleyOutlier}).
#' @export
#'
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Estimate the mean and covariance matrix
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
#' 
#' # Detect outliers using farness cutoff
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                      cutoff = "farness", 
#'                                      cutoff_lvl = 0.9)
#' 
#' # Compute Shapley values
#' credit_card_shapley <- int_Shapley(credit_card_int, 
#'                                    mean_c = credit_card_IMCD$mean_IMCD_c, 
#'                                    mean_r = credit_card_IMCD$mean_IMCD_r, 
#'                                    cov = credit_card_IMCD$cov_IMCD)
#' 
#' # Plot Shapley values with cutoff line and Interval-Mahalanobis distance
#' plot_bar_int_Shapley(credit_card_shapley, 
#'                     cutoff_value = credit_card_outliers$cutoff_value,
#'                     cutoff_label = "Farness 0.9",
#'                     palette = rainbow(credit_card_int@NIVar))
plot_bar_int_Shapley <- function(x, 
                                cutoff_value = NULL, 
                                cutoff_label = NULL, 
                                palette = NULL,
                                abbrev.var = 20, 
                                abbrev.obs = 20, 
                                sort.obs = TRUE, 
                                plot_IMah = TRUE,
                                IMah_label = expression(Robust~d[IMah]^2 * (bold(x))),
                                rotate_x = TRUE){

  if(!is.matrix(x)) stop("`x` must be a matrix.")

  phi <- as.data.frame(x)
  bar_order <- fill <- upper <- label <- x <- y <- NULL

  # Abbreviate names
  if(abbrev.obs > 0 && !is.null(rownames(phi))){
    rownames(phi) <- abbreviate(rownames(phi), minlength = abbrev.obs)
  }
  if(abbrev.var > 0 && !is.null(colnames(phi))){
    colnames(phi) <- abbreviate(colnames(phi), minlength = abbrev.var)
  }

  if(is.null(rownames(phi))){
    rownames(phi) <- paste0("Obs.", formatC(seq_len(nrow(phi)), width = nchar(nrow(phi)), flag = "0"))
  }
  if(is.null(colnames(phi))){
    colnames(phi) <- paste0("X", formatC(seq_len(ncol(phi)), width = nchar(ncol(phi)), flag = "0"))
  }

  rownames_vec <- rownames(phi)
  colnames_vec <- colnames(phi)

  # Sorting
  observation_sort <- if(sort.obs) {
    names(sort(rowSums(phi), decreasing = TRUE))
  } else rownames_vec

  # Palette
  if(is.null(palette)){
    if(!requireNamespace("RColorBrewer", quietly = TRUE)){
      stop("Package 'RColorBrewer' is required for palette functionality.")
    }
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    palette <- getPalette(max(12, ncol(phi)))
    palette <- palette[seq_len(ncol(phi))]
  }

  # Reshape
  plot_df <- as.data.frame(as.table(as.matrix(phi)))
  colnames(plot_df) <- c("rowname", "name", "value")

  plot_df$x <- factor(plot_df$rowname, levels = observation_sort)
  plot_df$fill <- factor(plot_df$name, levels = colnames_vec)
  plot_df$y <- plot_df$value

  plot_df$sign <- sign(plot_df$y)

  # Ordering within stacks
  plot_df <- plot_df[order(plot_df$x, -plot_df$sign, abs(plot_df$y)), ]
  plot_df$bar_order <- seq_len(nrow(plot_df))

  # IMah data
  if(plot_IMah){
    a_arrows <- data.frame(
      rowname = rownames_vec,
      lower = rowSums(phi * (phi < 0)),
      upper = rowSums(phi)
    )
    a_arrows$x <- factor(a_arrows$rowname, levels = observation_sort)
  }

  # Plot
  plt <- ggplot() +
    geom_bar(
      data = plot_df,
      aes(x = x, y = y, fill = fill, group = bar_order),
      stat = "identity", width = 0.9
    ) +
    geom_hline(yintercept = 0, linetype = 1) +
    labs(y = "Shapley Value") +
    theme_light() +
    theme(
      legend.box.margin = margin(-5, 0, -5, -40),
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "vertical",
      legend.direction = "horizontal",
      axis.title.x = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    ) +
    scale_fill_manual(name = NULL, values = palette)

  if(rotate_x){
    plt <- plt +
      theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
  }

  # IMah arrows
  if(plot_IMah){
    plt <- plt +
      geom_segment(
        data = a_arrows,
        aes(x = x, xend = x, y = 0, yend = upper, linetype = "IMah Distance"),
        linewidth = 1, color = "black",
        arrow = grid::arrow(angle = 90, length = grid::unit(0.1, "cm"),
                            ends = "last", type = "open")
      )
  }

  # Cutoff lines & labels (black lines, non-solid)
  if (!is.null(cutoff_value)) {
    if (!is.numeric(cutoff_value)) stop("cutoff_value must be numeric.")

    labels_vec <- if(!is.null(cutoff_label)) {
      if (length(cutoff_label) != length(cutoff_value)) {
        stop("cutoff_label must have same length as cutoff_value.")
      }
      cutoff_label
    } else if(!is.null(names(cutoff_value)) && any(names(cutoff_value) != "")) {
      names(cutoff_value)
    } else {
      paste0("Cutoff ", seq_along(cutoff_value))
    }

    linetypes <- rep(c("dashed","dotted","dotdash","longdash","twodash"),
                     length.out = length(cutoff_value))

    cutoff_df <- data.frame(
      y = cutoff_value,
      label = labels_vec,
      linetype = linetypes
    )

    plt <- plt +
      geom_hline(
        data = cutoff_df,
        aes(yintercept = y, linetype = label),
        color = "black", linewidth = 0.8
      ) +
      scale_linetype_manual(
        name = NULL,
        values = c(setNames(cutoff_df$linetype, cutoff_df$label),
                   "IMah Distance" = "solid"),
        labels = c(cutoff_df$label, IMah_label)
      )
  }

  plt
}

#' Tileplot of Shapley values for interval-valued data.
#'
#' @param shapley A \eqn{n \times p} matrix containing the Shapley values of \eqn{n} observations and \eqn{p} variables.
#' @param outliers A list containing the outliers' names as returned by \code{\link{int_outliers}}. If \code{outliers} is not \code{NULL} (default), only the outliers are highlighted in the plot.
#' @param rotate_x Logical. If \code{TRUE} (default), the x-axis labels are rotated.
#' @param abbrev.var Integer. If \code{abbrev.var} \eqn{> 0}, column names are abbreviated using abbreviate with \code{minlenght = abrev.var}.
#' @param abbrev.obs Integer. If \code{abbrev.obs} \eqn{> 0}, row names are abbreviated using abbreviate with \code{minlenght = abrev.obs}.
#' @param sort.var Logical. If \code{TRUE}, variables are sorted according to the distance.
#' @param sort.obs Logical. If \code{TRUE}, observations are sorted according to their squared Interval-Mahalanobis distance.
#' @param show_values Logical. If \code{TRUE}, the Shapley values are displayed in each tile.
#' @return Returns a tileplot that displays the Shapley values (\code{\link{int_Shapley}}) for each observation and variable. Optionally, only the outliers are highlighted in the plot.
#' @references Adapted from package \code{ShapleyOutlier} (\url{https://CRAN.R-project.org/package=ShapleyOutlier}).
#' @export
#' 
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Estimate the mean and covariance matrix
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
#' 
#' # Detect outliers using farness cutoff
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                      cutoff = "farness", 
#'                                      cutoff_lvl = 0.9)
#' 
#' # Compute Shapley values
#' credit_card_shapley <- int_Shapley(credit_card_int, 
#'                                    mean_c = credit_card_IMCD$mean_IMCD_c, 
#'                                    mean_r = credit_card_IMCD$mean_IMCD_r, 
#'                                    cov = credit_card_IMCD$cov_IMCD)
#' 
#' plot_tile_int_Shapley(credit_card_shapley, 
#'                      outliers = credit_card_outliers, 
#'                      sort.var = TRUE, 
#'                      sort.obs = TRUE)
plot_tile_int_Shapley <- function(shapley, 
                                 outliers = NULL, 
                                 rotate_x = TRUE, 
                                 abbrev.var = FALSE, 
                                 abbrev.obs = FALSE,
                                 sort.var = FALSE, 
                                 sort.obs = FALSE,
                                 show_values = FALSE) {

  if(!is.matrix(shapley)) stop("`shapley` must be a matrix.")

  Observation <- Feature <- Value <- Highlight <- NULL

  n <- nrow(shapley)
  p <- ncol(shapley)

  # Row/column names handling
  if (is.null(rownames(shapley))) {
    rownames(shapley) <- paste0(
      "Obs. ",
      formatC(seq_len(n), width = nchar(n), flag = "0")
    )
  }

  if (is.null(colnames(shapley))) {
    colnames(shapley) <- paste0(
      "X",
      formatC(seq_len(p), width = nchar(p), flag = "0")
    )
  }

  if (abbrev.obs > 0) {
    rownames(shapley) <- abbreviate(rownames(shapley), minlength = abbrev.obs)
  }

  if (abbrev.var > 0) {
    colnames(shapley) <- abbreviate(colnames(shapley), minlength = abbrev.var)
  }

  # Sorting
  if (sort.var) {
    features_sort <- colnames(shapley)[order(colSums(shapley))]
  } else {
    features_sort <- colnames(shapley)
  }

  if (sort.obs) {
    observation_sort <- rownames(shapley)[order(rowSums(shapley), decreasing = TRUE)]
  } else {
    observation_sort <- rownames(shapley)
  }

  # Base reshape
  df_melted <- as.data.frame(as.table(shapley))
  names(df_melted) <- c("Observation", "Feature", "Value")

  # Highlight handling
  if (is.null(outliers)) {
    df_melted$Highlight <- "Highlight"
  } else {
    df_melted$Highlight <- ifelse(
      df_melted$Observation %in% outliers$outliers_names,
      "Highlight",
      "NoHighlight"
    )
  }

  # Factor ordering
  df_melted$Observation <- factor(df_melted$Observation, levels = observation_sort)
  df_melted$Feature <- factor(df_melted$Feature, levels = features_sort)

  # Plot
  plt <- ggplot(df_melted, aes(
    x = Observation,
    y = Feature,
    fill = Value
  )) +
    geom_tile(
      color = "lightgray",
      aes(alpha = Highlight)
    ) +
    scale_fill_gradient2(
      low = "deepskyblue4",
      mid = "white",
      high = "red4",
      midpoint = 0
    ) +
    scale_alpha_manual(
      values = c("Highlight" = 1, "NoHighlight" = 0)
    ) +
    labs(fill = "Shapley\nValue") +
    theme_minimal() +
    guides(alpha = "none") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text = element_text(size = 22)
    )

  if (show_values) {
    plt <- plt +
      geom_text(
        aes(label = round(Value, 2)),
        size = 6
      ) +
      guides(fill = "none")
  }

  if (rotate_x) {
    plt <- plt +
      theme(
        axis.text.x = element_text(angle = 38, vjust = 1, hjust = 1)
      )
  }
  plt
}

#' Radar plot of Shapley values for interval-valued data.
#' 
#' @param shapley A \eqn{n \times p} matrix containing the Shapley values of \eqn{n} observations and \eqn{p} variables.
#' @param palette A vector of palette for each observation. Default is black.
#' @param sort.obs Logical. If \code{TRUE} (default), observations are sorted according to their squared (robust) Interval-Mahalanobis distance.
#' @return Returns a radar plot that displays the Shapley values (\code{\link{int_Shapley}}) for each observation.
#' @export
#'
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Estimate the mean and covariance matrix
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
#' 
#' # Detect outliers using farness cutoff
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                      cutoff = "farness", 
#'                                      cutoff_lvl = 0.9)
#' 
#' # Compute Shapley values
#' credit_card_shapley <- int_Shapley(credit_card_int, 
#'                                    mean_c = credit_card_IMCD$mean_IMCD_c,
#'                                    mean_r = credit_card_IMCD$mean_IMCD_r, 
#'                                    cov = credit_card_IMCD$cov_IMCD)
#' 
#' # colors
#' outliers_colors <- rep('black',credit_card_int@NObs)
#' names(outliers_colors) <- rownames(credit_card_int)
#' outliers_colors[credit_card_outliers$outliers_names] = '#009de0'
#' 
#' plot_radar_int_Shapley(credit_card_shapley, palette = outliers_colors)
plot_radar_int_Shapley <- function(shapley, 
                                  palette = NULL, 
                                  sort.obs = FALSE) {

  if(!is.matrix(shapley)) stop("`shapley` must be a matrix.")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  num_obs <- nrow(shapley)  # Number of observations

  if (is.null(palette)) {
    palette <- rep('black',num_obs)  # Default colors
  }

  # Order rows by rowSums(data)
  if(sort.obs){
    ordered_indices <- order(rowSums(shapley), decreasing = TRUE)
    shapley <- shapley[ordered_indices, , drop = FALSE]  # Reorder data
    palette <- palette[ordered_indices]  # Reorder colors accordingly
  }

  # Automatically determine rows & columns for the layout
  num_cols <- ceiling(sqrt(num_obs))  
  num_rows <- ceiling(num_obs / num_cols)  

  # Create min-max scaling reference
  radar_data <- rbind(
    max = apply(shapley, 2, max),
    min = apply(shapley, 2, min),
    as.data.frame(shapley)
  )

  # Define layout dynamically
  par(mfrow = c(num_rows, num_cols), mar = c(1, 1, 2, 1))

  # Generate radar plots for each observation
  for (i in 3:(num_obs + 2)) {
    fmsb::radarchart(
      df = radar_data[c(1, 2, i), ],
      pcol = palette[i - 2],  
      title = rownames(shapley)[i - 2],  
      pfcol = grDevices::adjustcolor(palette[i - 2], alpha.f = 0.5),  
      plwd = 2, plty = 1,
      cglcol = "grey", cglty = 1, cglwd = 0.8,
      axislabcol = "grey"
    )
  }
}

#' Beeswarm plot of Shapley values for interval-valued data.
#' 
#' @param shapley A \eqn{n \times p} matrix containing the Shapley values of \eqn{n} observations and \eqn{p} variables.
#' @param color_class A vector indicating the color class of each observation. If NULL (default), all points have the same color.
#' @param color_label Character. Label for the color class. If NULL (default), no legend for the color class is shown.
#' @param palette A vector with colors for each color class. Default is NULL.
#' @param rotate_x Logical. If \code{TRUE} (default), the x-axis labels are rotated.
#' @param shape_class A vector indicating the shape class of each observation. If NULL (default), all points have the same shape.
#' @param shape_label Character. Label for the shape class. If NULL (default), no legend for the shape class is shown.
#' @param ggplotly Logical. If \code{TRUE} (default), the plot is converted to an interactive [plotly] object.
#' @param label_obs A vector with the names of the observations to be labeled in the plot when \code{ggplotly = FALSE}. Default is NULL.
#' @return Returns a beeswarm plot that displays the Shapley values (\code{\link{int_Shapley}}) for each observation and feature.
#' @export
#' 
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Estimate the mean and covariance matrix
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
#' 
#' # Detect outliers using farness cutoff
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                      cutoff = "farness", 
#'                                      cutoff_lvl = 0.9)
#' 
#' # Compute Shapley values
#' credit_card_shapley <- int_Shapley(credit_card_int, 
#'                                    mean_c = credit_card_IMCD$mean_IMCD_c,
#'                                    mean_r = credit_card_IMCD$mean_IMCD_r, 
#'                                    cov = credit_card_IMCD$cov_IMCD)
#' 
#' # Beeswarm plot of Shapley values colored by outlier status
#' plot_beeswarm_int_Shapley(credit_card_shapley, 
#'                      color_class = credit_card_outliers$is_outlier, 
#'                      palette = c("gray50", "darkred"), 
#'                      color_label = "Outlier Status")
plot_beeswarm_int_Shapley <- function(shapley, 
                                 color_class, 
                                 color_label = NULL, 
                                 palette = NULL, 
                                 rotate_x = TRUE,
                                 shape_class = NULL, 
                                 shape_label = NULL, 
                                 ggplotly = FALSE, 
                                 label_obs = NULL) {

  if(!is.matrix(shapley)) stop("`shapley` must be a matrix.")

  Feature <- Shapley_Value <- Color_Class <- Shape_Class <- Observation <- NULL

  # Convert to long format
  shap_long <- as.data.frame(as.table(shapley))
  names(shap_long) <- c("Observation", "Feature", "Shapley_Value")

  # Ensure classes are valid
  if (length(color_class) != nrow(shapley)) {
    stop("`color_class` must have length equal to number of rows of `shapley`.")
  }

  # Base data with color class
  base_data <- data.frame(
    Observation = rownames(shapley),
    Color_Class = color_class,
    stringsAsFactors = FALSE
  )

  # Add shape class if provided
  if (!is.null(shape_class)) {
    if (length(shape_class) != nrow(shapley)) {
      stop("`shape_class` must have length equal to number of rows of `shapley`.")
    }
    base_data$Shape_Class <- shape_class
  }

  # Merge shapley data with class info
  shap_data <- cbind(shap_long, base_data[match(shap_long$Observation, base_data$Observation), ])

  # Feature ordering
  shap_data$Feature <- factor(shap_data$Feature, levels = colnames(shapley))

  # Set up ggplot aesthetics depending on shape_class
  if (!is.null(shape_class)) {
    p <- ggplot(
      shap_data,
      aes(
        x = Feature,
        y = Shapley_Value,
        color = Color_Class,
        shape = Shape_Class,
        text = Observation
      )
    )
  } else {
    p <- ggplot(
      shap_data,
      aes(
        x = Feature,
        y = Shapley_Value,
        color = Color_Class,
        text = Observation
      )
    )
  }

  # Base beeswarm plot
  p <- p +
    ggbeeswarm::geom_quasirandom(
      data = subset(shap_data, Color_Class == "Regular"),
      width = 0.2,
      size = 1.5,
      alpha = 0.7,
      varwidth = TRUE
    ) +
    ggbeeswarm::geom_quasirandom(
      data = subset(shap_data, Color_Class != "Regular"),
      width = 0.2,
      size = 1.5,
      alpha = 0.7,
      varwidth = TRUE
    ) +
    theme_bw() +
    labs(
      y = "Shapley Value",
      color = color_label,
      shape = shape_label
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")
    )

  # Palette
  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }

  # Conditionally hide legends
  if (is.null(color_label) && (is.null(shape_label) || is.null(shape_class))) {
    p <- p + theme(legend.position = "none")
  } else {
    if (is.null(color_label)) {
      p <- p + guides(color = "none") +
        theme(legend.position = "top", 
                        legend.text = element_text(size = 12),
                        legend.box.margin = margin(-5, 0, -5, -10),
                        legend.margin = margin(0, 0, 0, 0))
    }
    if (is.null(shape_label) || is.null(shape_class)) {
      p <- p + guides(shape = "none") +
        theme(legend.position = "top", 
                        legend.text = element_text(size = 12),
                        legend.box.margin = margin(-5, 0, -5, 0),
                        legend.margin = margin(0, 0, 0, 0))
    }
  }

  # Add vertical separator lines between features
  feature_levels <- colnames(shapley)

  if (length(feature_levels) > 1) {
    vline_positions <- seq_along(feature_levels)[-length(feature_levels)] + 0.5

    p <- p + geom_vline(
      xintercept = vline_positions,
      color = "grey70"
    )
  }

  # Optional labeling of specific observations
  if (!ggplotly && !is.null(label_obs)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("Package 'ggrepel' is required for labeling functionality.")
    }
    df_labeled <- subset(shap_data, Observation %in% label_obs)
    if (nrow(df_labeled) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = df_labeled,
          aes(label = Observation, color = Color_Class),
          size = 4,
          max.overlaps = Inf,
          box.padding = 0.3,
          point.padding = 0.3
        )
    }
  }

  # Rotate x labels
  if (rotate_x) {
    p <- p +
      theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
  }

  # Output
  if (ggplotly) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive output.")
    }
    plotly::ggplotly(p, tooltip = c("y", "text"))
  } else {
    p
  }
}

#' Barplot of Shapley value decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) for interval-valued data.
#' 
#' @param shapley_decomp A list of matrices containing the Shapley value decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) for each observation.
#' @param palette A vector with colors for each feature. If \code{palette} is \code{NULL} (default), the colors are generated using \code{RColorBrewer}.
#' @param rotate_x Logical. If \code{TRUE} (default), the x-axis labels are rotated.
#' @param abbrev.obs Integer. If \code{abbrev.obs} \eqn{> 0}, row names are abbreviated using abbreviate with \code{minlenght = abbrev.obs}.
#' @param sort.obs Logical. If \code{TRUE} (default), observations are sorted according to their total Shapley value.
#' @param plot_IMah Logical. If \code{TRUE}, the Interval-Mahalanobis distance (sum of all Shapley values) will be included in the plot.
#' @return Returns a barplot that displays the Shapley value decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) for each observation.
#' @export
#' 
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute Shapley decomposition into contributions of Centers, Ranges, and CrossCentersRanges
#' # based on IMCD estimates of mean and covariance matrix
#' credit_card_shap_decomp <- int_Shapley_decomp(credit_card_int)
#'
#' # Plot Shapley decomposition with contributions of Centers, Ranges, and CrossCentersRanges
#' plot_bar_int_Shapley_decomp(credit_card_shap_decomp, palette = rainbow(credit_card_int@NIVar))
plot_bar_int_Shapley_decomp <- function(shapley_decomp, 
                                       palette = NULL, 
                                       rotate_x = TRUE, 
                                       abbrev.obs = 20,
                                       sort.obs = TRUE,
                                       plot_IMah = FALSE) {

  if (!is.list(shapley_decomp) || !all(vapply(shapley_decomp, is.matrix, logical(1)))) {
    stop("`shapley_decomp` must be a list of matrices.")
  }

  Row <- Column <- Value <- y <- label <- NULL  # Avoid warnings
  
  # Abbreviate names
  if (abbrev.obs > 0 && !is.null(names(shapley_decomp))) {
    names(shapley_decomp) <- abbreviate(names(shapley_decomp), minlength = abbrev.obs)
  }

  # Sort observations
  if (sort.obs) {
    ord <- order(vapply(shapley_decomp, sum, numeric(1)), decreasing = TRUE)
    shapley_decomp <- shapley_decomp[ord]
  }

  # Palette
  if(is.null(palette)){
    if(!requireNamespace("RColorBrewer", quietly = TRUE)){
      stop("Package 'RColorBrewer' is required for palette functionality.")
    }
    colorCount <- ncol(shapley_decomp[[1]])
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    palette <- if (colorCount > 12) getPalette(colorCount) else getPalette(12)
  }

  mat_names <- names(shapley_decomp)
  if (is.null(mat_names)) mat_names <- seq_along(shapley_decomp)

  # Convert list of matrices to long format data frame
  df_long <- do.call(rbind, lapply(seq_along(shapley_decomp), function(i) {
    m <- shapley_decomp[[i]]
    rn <- rownames(m)
    cn <- colnames(m)

    transform(
      expand.grid(Row = rn,
                  Column = cn,
                  KEEP.OUT.ATTRS = FALSE,
                  stringsAsFactors = FALSE),
      Value = as.vector(m),
      Matrix = mat_names[i]
    )
  }))

  df_long$Column <- factor(df_long$Column, levels = colnames(shapley_decomp[[1]]))
  df_long$Row <- factor(df_long$Row, levels = rownames(shapley_decomp[[1]]))
  df_long$Matrix <- factor(df_long$Matrix, levels = mat_names)

  plt <- ggplot(df_long, aes(x = Row, y = Value, fill = Column)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Matrix, ncol = length(shapley_decomp)) +
    labs(y = "Shapley Value", fill = "Features") +
    theme_light() +
    theme(
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(size = 12, face = "bold", colour = "black"),
      legend.box.margin = margin(-5, 0, -5, -10),
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = "top",
      legend.justification = "center",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12)
    ) +
    scale_fill_manual(values = palette) +
    guides(fill = guide_legend(ncol = length(shapley_decomp)-2, byrow = TRUE))

  if (plot_IMah) {
    IMah_df <- data.frame(
      Matrix = factor(mat_names, levels = mat_names),
      y = vapply(shapley_decomp, function(m) sum(m), numeric(1)),
      label = "Robust Squared Distance"
    )

    plt <- plt +
      geom_hline(
        data = IMah_df,
        aes(yintercept = y, linetype = label),
        linewidth = 0.8,
        colour = "black"
      ) +
      scale_linetype_manual(
        values = c("Robust Squared Distance" = "dashed"),
        name = NULL
      )
  }

  if (rotate_x) {
    plt <- plt +
      theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
  }

  plt
}
