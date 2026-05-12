#' Symbolic Biplot for Interval-valued Data
#' 
#' Create a biplot for interval-valued symbolic data, visualizing the symbolic data as rectangles or crosses, with the first two variables on the x and y axes. The function allows customization of colors, fill colors, and outlier representation.
#' 
#' @param data An \linkS4class{intData} object containing the macrodata/interval data. The first two variables are used for the x and y axes.
#' @param type The type of plot to generate: "rectangles", "crosses" or "crosses2". Default is "rectangles".
#' @param palette A vector with colors for each observation. Default is \code{rainbow(nrow(data))}.
#' @param fill_col If \code{type="rectangles"}, a vector with colors for the fill of each observation, or a single color for all observations. Default is "gray50".
#' @param is_outlier A vector with logical values indicating if the observation is an outlier or not. It makes the line width of the outlying observations thicker. Default is NULL.
#' @param ... Additional graphical parameters.
#' @return A biplot is drawn in the graphic window. The biplot shows the symbolic data as rectangles or crosses, with the first two variables on the x and y axes.
#' @importFrom graphics rect lines
#' @importFrom grDevices adjustcolor rainbow
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' SYMB.biplot(credit_card_int[,c(3,5)])
#' 
#' # Highlight outliers in the biplot
#' credit_card_IMCD <- IMCD(credit_card_int, floor(0.75*credit_card_int@NObs), "farness", 0.9)
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, "farness", 0.9)
#' outliers_colors<-rep('gray50',credit_card_int@NObs)
#' names(outliers_colors)<-rownames(credit_card_int)
#' outliers_colors[credit_card_outliers$outliers_names] = 'red'
#' SYMB.biplot(credit_card_int[,c(3,5)], palette = outliers_colors, 
#'             is_outlier = credit_card_outliers$is_outlier)
#' @export
SYMB.biplot <- function(data,
                        type = c("rectangles", "crosses", "crosses2"),
                        palette = rainbow(nrow(data)),
                        fill_col = "gray50",
                        is_outlier = NULL,
                        ...) {
  type <- match.arg(type)

  Lbnd <- LowerBounds(data)
  Ubnd <- UpperBounds(data)

  xmin <- Lbnd[, 1]; xmax <- Ubnd[, 1]
  ymin <- Lbnd[, 2]; ymax <- Ubnd[, 2]

  if (length(fill_col) == 1) fill_col <- rep(fill_col, data@NObs)
  if (is.null(is_outlier)) is_outlier <- rep(FALSE, data@NObs)

  # Create base plot
  plot(range(c(xmin, xmax)), range(c(ymin, ymax)), 
        type = "n", xlab = colnames(data)[1], ylab = colnames(data)[2], ...)

  draw_symbol <- function(i, lwd = 1) {
    col_border <- palette[i]
    col_fill <- adjustcolor(fill_col[i], alpha.f = 0.3)

    if (type == "rectangles") {
      rect(xmin[i], ymin[i], xmax[i], ymax[i],
           border = col_border, col = col_fill, lwd = lwd, ...)
    }

    if (type == "crosses2") {
      lines(c(xmin[i], xmax[i]), c(ymin[i], ymax[i]), col = col_border, lwd = lwd, ...)
      lines(c(xmin[i], xmax[i]), c(ymax[i], ymin[i]), col = col_border, lwd = lwd, ...)
    }

    if (type == "crosses") {
      ymid <- (ymin[i] + ymax[i]) / 2
      xmid <- (xmin[i] + xmax[i]) / 2
      lines(c(xmin[i], xmax[i]), rep(ymid, 2), col = col_border, lwd = lwd, ...)
      lines(rep(xmid, 2), c(ymin[i], ymax[i]), col = col_border, lwd = lwd, ...)
    }
  }

  # Draw non-outliers first
  for (i in which(!is_outlier)) {
    draw_symbol(i, lwd = 1)
  }

  # Draw outliers on top
  for (i in which(is_outlier)) {
    draw_symbol(i, lwd = 2)
  }
}


#' Pairs-plot for Interval-valued Symbolic data.
#'
#' Adapted from pairs.panels (R package "psych") shows a scatter plot of matrices, with bivariate symbolic scatter plots below the diagonal, variables' names on the diagonal, and all the symbolic correlations above the diagonal. Useful for descriptive statistics of symbolic objects described by interval variables.
#' 
#' @param data An \linkS4class{intData} object containing the macrodata/interval data
#' @param type The type of plot to generate: "rectangles" or "crosses" or "crosses2". Default is "rectangles".
#' @param cex.cor Character expansion factor
#' @param corr A matrix with the symbolic correlations; if not provided the upper panel is omitted
#' @param palette A vector with colors for each observation.
#' @param fill_col If \code{type="rectangles"}, a vector with colors for the fill of each observation, or a single color for all observations. Default is "gray50".
#' @param is_outlier A vector with logical values indicating if the observation is an outlier or not. It makes the line width of the outlying observations thicker. Default is NULL.
#' @param ... Additional graphical parameters.
#' @return A scatter plot matrix is drawn in the graphic window. The lower off diagonal draws scatter plots, the diagonal variables' names, the upper off diagonal reports  all the symbolic correlations.
#' @importFrom graphics pairs rect lines text par
#' @importFrom grDevices adjustcolor colorRampPalette rainbow
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_cov<-int_cov(credit_card_int)
#' credit_card_cor<-cov2cor(credit_card_cov)
#' SYMB.pairs.panels(credit_card_int,corr=credit_card_cor,labels=colnames(credit_card_int))
#' 
#' # Highlight outliers in the biplot
#' credit_card_IMCD <- IMCD(credit_card_int, floor(0.75*credit_card_int@NObs), "farness", 0.9)
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, "farness", 0.9)
#' outliers_colors<-rep('gray50',credit_card_int@NObs)
#' names(outliers_colors)<-rownames(credit_card_int)
#' outliers_colors[credit_card_outliers$outliers_names] = 'red'
#' SYMB.pairs.panels(credit_card_int,corr=cov2cor(credit_card_IMCD$cov_IMCD), 
#'                  palette = outliers_colors,labels=colnames(credit_card_int),
#'                  type = "rectangles",is_outlier = credit_card_outliers$is_outlier)
#' @export
SYMB.pairs.panels<-function (data,
                              type=c("rectangles","crosses","crosses2"),
                              cex.cor=2.0,
                              corr=NULL,
                              palette=rainbow(nrow(data)),
                              fill_col="gray50",
                              is_outlier=NULL,
                              ...){
  type <- match.arg(type)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  "SYMB.panel.rect" <- function(x,y,n=length(x)/2, ...) {
      xmin<-x[1:(length(x)/2)]
      xmax<-x[(length(x)/2+1):length(x)]
      ymin<-y[1:(length(y)/2)]
      ymax<-y[(length(y)/2+1):length(y)]

      if (length(fill_col) == 1) fill_col <- rep(fill_col, n)
      if (is.null(is_outlier)) is_outlier <- rep(FALSE, n)

      draw_symbol <- function(i, lwd = 1) {
        col_border <- palette[i]
        col_fill <- adjustcolor(fill_col[i], alpha.f = 0.3)

        if (type == "rectangles") {
          rect(xmin[i], ymin[i], xmax[i], ymax[i],
              border = col_border, col = col_fill, lwd = lwd, ...)
        }

        if (type == "crosses2") {
          lines(c(xmin[i], xmax[i]), c(ymin[i], ymax[i]), col = col_border, lwd = lwd, ...)
          lines(c(xmin[i], xmax[i]), c(ymax[i], ymin[i]), col = col_border, lwd = lwd, ...)
        }

        if (type == "crosses") {
          ymid <- (ymin[i] + ymax[i]) / 2
          xmid <- (xmin[i] + xmax[i]) / 2
          lines(c(xmin[i], xmax[i]), rep(ymid, 2), col = col_border, lwd = lwd, ...)
          lines(rep(xmid, 2), c(ymin[i], ymax[i]), col = col_border, lwd = lwd, ...)
        }
      }

      # Draw non-outliers first
      for (i in which(!is_outlier)) {
        draw_symbol(i, lwd = 1)
      }

      # Draw outliers on top
      for (i in which(is_outlier)) {
        draw_symbol(i, lwd = 2)
      }
    }

  "SYMB.panel.cor" <- function(x, y, digits = 3, ...){
      par(usr = c(0, 1, 0, 1))
      i <- which(colSums(SYMB.matrix==x)==nrow(SYMB.matrix))
      j <- which(colSums(SYMB.matrix==y)==nrow(SYMB.matrix))
      corr_value <- corr[i,j]
      txt <- format(round(corr_value, digits), nsmall = digits)
      color_palette <- colorRampPalette(c("deepskyblue4", "white", "red4"))
      color <- color_palette(100)[findInterval(corr_value, seq(-1, 1, length.out = 100))]
      rect(0, 0, 1, 1, col = color)
      text(0.5, 0.5, txt, cex = cex.cor)
  }
  Lbnd<-LowerBounds(data); Ubnd<-UpperBounds(data)
  names(Lbnd)<-names(Ubnd)<-colnames(data)
  SYMB.matrix<-rbind(Lbnd,Ubnd)
  if (is.null(corr)){
      pairs(SYMB.matrix, upper.panel = NULL, lower.panel = SYMB.panel.rect,  ...)
  }else{
      pairs(SYMB.matrix, upper.panel = SYMB.panel.cor, lower.panel = SYMB.panel.rect,  ...)
  }
}

#' Distance-Distance plot for interval-valued data.
#'
#' @param class_dist A numeric vector containing the classical distances for each observation.
#' @param class_cutoff Numeric. The cutoff value for the classical distances.
#' @param class_cutoff_label Character. Label for the classical cutoff. If NULL (default), no legend for the classical cutoff is shown.
#' @param rob_dist A numeric vector containing the robust distances for each observation.
#' @param rob_cutoff Numeric. The cutoff value for the robust distances.
#' @param rob_cutoff_label Character. Label for the robust cutoff. If NULL (default), no legend for the robust cutoff is shown.
#' @param obs_names A character vector containing the names of the observations. If NULL (default), the names are taken from the names of class_dist.
#' @param ggplotly Logical. If \code{TRUE} (default), the plot is converted to an interactive [plotly] object.
#' @param color_class A vector indicating the color class of each observation. If NULL (default), all points have the same color.
#' @param color_label Character. Label for the color class. If NULL (default), no legend for the color class is shown.
#' @param palette A vector with colors for each color class. If NULL (default), default [ggplot2] colors are used.
#' @param shape_class A vector indicating the shape class of each observation. If NULL (default), all points have the same shape.
#' @param shape_label Character. Label for the shape class. If NULL (default), no legend for the shape class is shown.
#' @param label_obs A vector with the names of the observations to be labeled in the plot when \code{ggplotly = FALSE}. Default is NULL.
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom ggrepel geom_text_repel
#' @return Returns a Distance-Distance plot that displays the classical distances against the robust distances for each observation, highlighting outliers.
#' @export
#' @examples
#' #Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' #Estimate the mean and covariance matrix
#' credit_card_IMCD<-IMCD(credit_card_int, floor(nrow(credit_card_int)*0.75), "farness", 0.9)
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                            p=credit_card_int@NIVar, cutoff_lvl = 0.9)
#' 
#' #Plot Distance-Distance plot
#' class_dist <- IMah_dist(credit_card_int, z=rep(1,credit_card_int@NObs))
#' class_outliers <- int_outliers(class_dist,cutoff = "adjbox",p=p,cutoff_lvl = 1.5)
#' credit_card_is_outliers <- as.character(credit_card_outliers$is_outlier)
#' credit_card_is_outliers[credit_card_outliers$is_outlier] <- "Outlier"
#' credit_card_is_outliers[!credit_card_outliers$is_outlier] <- "Inlier"
#' plot_dist_dist(class_dist, class_outliers$cutoff_value[2], "1.5 adjusted boxplot",
#'               credit_card_IMCD$robust_dist, credit_card_outliers$cutoff_value, "0.9 farness",
#'               color_class = credit_card_is_outliers, palette = c("grey50", "red"))
plot_dist_dist <- function(class_dist, 
                          class_cutoff = NULL, 
                          class_cutoff_label = NULL,
                          rob_dist, 
                          rob_cutoff = NULL, 
                          rob_cutoff_label = NULL,
                          obs_names = NULL, 
                          ggplotly = TRUE,
                          color_class = NULL, 
                          color_label = NULL, 
                          palette = NULL,
                          shape_class = NULL, 
                          shape_label = NULL,
                          label_obs = NULL){
  Classical <- Robust <- Name <- Color_Class <- Shape_Class  <- x  <- y  <- type <- NULL

  if (is.null(obs_names)) {
    obs_names <- names(class_dist)
  }

  df <- data.frame(
    Classical = class_dist,
    Robust    = rob_dist,
    Name      = obs_names
  )

  # Color class
  if (!is.null(color_class)) {
    if (length(color_class) != nrow(df))
      stop("`color_class` must match length of class_dist.")
    df$Color_Class <- as.factor(color_class)
  }

  # Shape class
  if (!is.null(shape_class)) {
    if (length(shape_class) != nrow(df))
      stop("`shape_class` must match length of class_dist.")
    df$Shape_Class <- as.factor(shape_class)
  }

  p <- ggplot(df, aes(x = Classical, y = Robust, label = Name))

  if (!is.null(color_class) && !is.null(shape_class)) {
    p <- p + geom_point(aes(color = Color_Class, shape = Shape_Class), size = 3)
  } else if (!is.null(color_class)) {
    p <- p + geom_point(aes(color = Color_Class), size = 3)
  } else if (!is.null(shape_class)) {
    p <- p + geom_point(aes(shape = Shape_Class), size = 3)
  } else {
    p <- p + geom_point(size = 3)
  }

  # Provide default labels if missing
  if (!is.null(class_cutoff) && is.null(class_cutoff_label))
    class_cutoff_label <- paste("Classical cutoff", seq_along(class_cutoff))

  if (!is.null(rob_cutoff) && is.null(rob_cutoff_label))
    rob_cutoff_label <- paste("Robust cutoff", seq_along(rob_cutoff))

  # Create DF for cutoffs
  cutoff_df <- data.frame()

  if (!is.null(class_cutoff)) {
    cutoff_df <- rbind(
      cutoff_df,
      data.frame(
        x = class_cutoff,
        y = NA,
        type = factor(class_cutoff_label, levels = c(class_cutoff_label, rob_cutoff_label)),
        axis = "x"
      )
    )
  }

  if (!is.null(rob_cutoff)) {
    cutoff_df <- rbind(
      cutoff_df,
      data.frame(
        x = NA,
        y = rob_cutoff,
        type = factor(rob_cutoff_label, levels = c(class_cutoff_label, rob_cutoff_label)),
        axis = "y"
      )
    )
  }

  # Assign linetypes automatically
  if (nrow(cutoff_df) > 0) {
    unique_types <- unique(cutoff_df$type)
    auto_linetypes <- rep(c("dashed","dotted","dotdash","twodash","longdash","F1"), length.out = length(unique_types))
    names(auto_linetypes) <- unique_types

    # Add classical cutoffs
    if (!is.null(class_cutoff)) {
      p <- p +
        geom_vline(
          data = cutoff_df[cutoff_df$axis == "x", ],
          aes(xintercept = x, linetype = type),
          color = "black"
        )
    }

    # Add robust cutoffs
    if (!is.null(rob_cutoff)) {
      p <- p +
        geom_hline(
          data = cutoff_df[cutoff_df$axis == "y", ],
          aes(yintercept = y, linetype = type),
          color = "black"
        )
    }

    # Add linetype guide
    p <- p + scale_linetype_manual(name = "Cutoff", values = auto_linetypes)
  }

  # Labels & Theme
  p <- p +
    labs(
      x = "Squared Classical Distance",
      y = "Squared Robust Distance",
      color = color_label,
      shape = shape_label
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.box = "vertical",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      legend.box.margin = margin(-5, 0, -5, 0),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 14, face = "bold")
    )

  if (!is.null(palette) && !is.null(color_class)) {
    p <- p + scale_color_manual(values = palette)
  }

  # Conditional legend hiding
  if (is.null(color_label)) p <- p + guides(color = "none")
  if (is.null(shape_label) || is.null(shape_class)) p <- p + guides(shape = "none")

  # Optional labeling of specific observations
  if (!ggplotly && !is.null(label_obs)) {
    df_labeled <- subset(df, Name %in% label_obs)
    if (nrow(df_labeled) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = df_labeled,
          aes(label = Name, color = Color_Class),
          size = 4,
          max.overlaps = Inf,
          box.padding = 0.3,
          point.padding = 0.3
        )
    }
  }

  # Convert to plotly if requested
  if (ggplotly) {
    plotly::ggplotly(p, tooltip = c("label", "x", "y"))
  } else {
    p
  }
}


#' Interval-Mahalanobis distance plot for interval-valued data.
#'
#' @param dist A numeric vector containing the Interval-Mahalanobis distances for each observation.
#' @param cutoff A numeric vector containing cutoff values to be displayed as horizontal lines.
#' @param cutoff_label A character vector containing labels for each cutoff. If NULL (default), default labels are generated.
#' @param obs_names A character vector containing the names of the observations. If NULL (default), the names are taken from the names of dist.
#' @param sort.obs Logical. If \code{TRUE} (default), observations are sorted according to their distances.
#' @param color_class A vector indicating the color class of each observation. If NULL (default), all points have the same color.
#' @param color_label Character. Label for the color class. If NULL (default), no legend for the color class is shown.
#' @param palette A vector with colors for each color class. If NULL (default), default [ggplot2] colors are used.
#' @param shape_class A vector indicating the shape class of each observation. If NULL (default), all points have the same shape.
#' @param shape_label Character. Label for the shape class. If NULL (default), no legend for the shape class is shown.
#' @param label_obs A vector with the names of the observations to be labeled in the plot. If NULL (default), no labels are shown and x-axis labels are displayed.
#' @return Returns a plot that displays the Interval-Mahalanobis distances for each observation, highlighting outliers based on specified cutoffs.
#' @export
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats setNames
#' @examples
#' #Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' #Estimate the mean and covariance matrix
#' credit_card_IMCD<-IMCD(credit_card_int, floor(nrow(credit_card_int)*0.75), "farness", 0.9)
#' credit_card_outliers <- int_outliers(credit_card_IMCD$robust_dist, 
#'                                            p=credit_card_int@NIVar, cutoff_lvl = 0.9)
#' credit_card_is_outliers <- as.character(credit_card_outliers$is_outlier)
#' credit_card_is_outliers[credit_card_outliers$is_outlier] <- "Outlier"
#' credit_card_is_outliers[!credit_card_outliers$is_outlier] <- "Inlier"
#' 
#' #Plot Interval-Mahalanobis distance plot
#' plot_interval_dist(credit_card_IMCD$robust_dist,
#'                    cutoff = credit_card_outliers$cutoff_value,
#'                    cutoff_label = c("0.9 farness"),
#'                    obs_names = rownames(credit_card_int),
#'                    sort.obs = FALSE,
#'                    color_class = credit_card_is_outliers,
#'                    palette = c("grey50", "red"))
plot_interval_dist <- function(
  dist,
  cutoff = NULL,
  cutoff_label = NULL,
  obs_names = NULL,
  sort.obs = TRUE,
  color_class = NULL,
  color_label = NULL,
  palette = NULL,
  shape_class = NULL,
  shape_label = NULL,
  label_obs = NULL
) {

  obs <- value <- Color_Class <- Shape_Class <- type <- y <- NULL

  if (is.null(obs_names)) {
    obs_names <- names(dist)
  }

  df <- data.frame(
    obs   = obs_names,
    value = dist,
    stringsAsFactors = FALSE
  )

  # Sort observations (and aesthetics!)
  if (sort.obs) {
    ord <- order(df$value, decreasing = TRUE)
    df <- df[ord, ]

    if (!is.null(color_class)) color_class <- color_class[ord]
    if (!is.null(shape_class)) shape_class <- shape_class[ord]
  }

  df$obs <- factor(df$obs, levels = df$obs)

  # Color class
  if (!is.null(color_class)) {
    if (length(color_class) != nrow(df))
      stop("`color_class` must match length of dist.")
    df$Color_Class <- as.factor(color_class)
  }

  # Shape class
  if (!is.null(shape_class)) {
    if (length(shape_class) != nrow(df))
      stop("`shape_class` must match length of dist.")
    df$Shape_Class <- as.factor(shape_class)
  }

  p <- ggplot(df, aes(x = obs, y = value))

  if (!is.null(color_class) && !is.null(shape_class)) {
    p <- p + geom_point(aes(color = Color_Class, shape = Shape_Class), size = 3)
  } else if (!is.null(color_class)) {
    p <- p + geom_point(aes(color = Color_Class), size = 3)
  } else if (!is.null(shape_class)) {
    p <- p + geom_point(aes(shape = Shape_Class), size = 3)
  } else {
    p <- p + geom_point(size = 3)
  }

  # Cutoffs
  if (!is.null(cutoff)) {

    if (is.null(cutoff_label)) {
      cutoff_label <- paste(
        "Interval-valued Mahalanobis Distance cutoff",
        seq_along(cutoff)
      )
    }

    cutoff_df <- data.frame(
      y = cutoff,
      type = factor(cutoff_label, levels = cutoff_label)
    )

    auto_linetypes <- rep(
      c("dashed","dotted","dotdash","twodash","longdash"),
      length.out = length(cutoff_label)
    )
    names(auto_linetypes) <- cutoff_label

    p <- p +
      geom_hline(
        data = cutoff_df,
        aes(yintercept = y, linetype = type),
        color = "black",
        linewidth = 1
      ) +
      scale_linetype_manual(name = "Cutoff", values = auto_linetypes)
  }

  # Apply palette safely
  if (!is.null(palette) && !is.null(color_class)) {
    levs <- levels(df$Color_Class)
    if (is.null(names(palette)) || !all(levs %in% names(palette))) {
      palette <- setNames(palette[seq_along(levs)], levs)
    }
    p <- p + scale_color_manual(values = palette)
  }

  # Optional ggrepel labels
  if (!is.null(label_obs)) {
    label_df <- df[df$obs %in% label_obs, , drop = FALSE]

    if (nrow(label_df) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = obs, color = Color_Class),
          size = 4,
          max.overlaps = Inf,
          box.padding = 0.3,
          point.padding = 0.3,
          show.legend = FALSE
        )
    }
  }

  # Base theme
  p <- p +
    labs(
      x = NULL,
      y = "Squared Robust Distance",
      color = color_label,
      shape = shape_label
    ) +
    theme_light() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "top",
      legend.box = "vertical",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      legend.box.margin = margin(-5, 0, -5, 0),
      legend.margin = margin(0, 0, 0, 0)
    )

  # CONDITIONAL axis behavior
  if (is.null(label_obs)) {
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 25)
    )
  } else {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  if (is.null(color_label)) p <- p + guides(color = "none")
  if (is.null(shape_label) || is.null(shape_class)) p <- p + guides(shape = "none")

  p
}
