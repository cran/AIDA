#' Farness Estimation
#' 
#' Estimate farness from a distance vector in order to identify outlier observations.
#' 
#' @param dist Vector of distances of each observation.
#' @param cutoff_value Optional cutoff value between 0 and 1 to flag outliers. If provided, the function returns both the farness probabilities and the cutoff distance value in the original distance scale.
#' @return  Farness of each observation. Values between 0 and 1. If \code{cutoff_value} is provided, a list with the farness probabilities and the cutoff distance value in the original distance scale is returned.
#' @importFrom cellWise transfo transfo_transformback
#' @importFrom stats mad pnorm qnorm sd median
#' @references J. Raymaekers and P.J. Rousseeuw (2021). Transforming variables to central normality. Machine Learning. \doi{10.1007/s10994-021-05960-5}
#' @references Based on the \code{cellWise} package: Raymaekers J, Rousseeuw P (2023). _cellWise: Analyzing Data with Cellwise Outliers_. R package version 2.5.3, \url{https://CRAN.R-project.org/package=cellWise}.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute squared Interval-Mahalanobis distance
#' z <- rep(1, nrow(credit_card_int))
#' credit_card_dist<-IMah_dist(credit_card_int,z)
#' 
#' credit_card_farness <- farness(credit_card_dist, 0.9)
farness <- function(dist, cutoff_value=NULL) {
    indnz <- which(dist > 1e-10)
    farnz <- dist[indnz]

    #scale the distances
    farloc <- median(farnz, na.rm = TRUE)
    farsca <- mad(farnz, na.rm = TRUE)
    if (farsca < 1e-10) {farsca <- sd(farnz, na.rm = TRUE)}
    sfar <- scale(farnz, center = farloc, scale = farsca)

    #apply Yeo-Johnson transformation
    YJout  <- cellWise::transfo(X = sfar, type = "YJ", robust = TRUE, standardize = FALSE, checkPars = list(silent = TRUE))
    xt <- YJout$Y

    #scale the transformed distances
    tfarloc <- median(xt, na.rm = TRUE)
    tfarsca <- mad(xt, na.rm = TRUE)
    zt <- scale(xt, tfarloc, tfarsca)
    
    #obtain final probabilities
    probs <- rep(0, length(dist))
    probs[indnz] <- pnorm(zt)
    
    #cutoff value
    if (!is.null(cutoff_value)) {
        if (cutoff_value>1||cutoff_value<0) stop("cutoff_value must be between 0 and 1")
        cutoff_value <- qnorm(cutoff_value)
        cutoff_value <- cutoff_value * tfarsca + tfarloc
        cutoff_value <- cellWise::transfo_transformback(Ynew = cutoff_value, transfo.out = YJout)
        cutoff_value <- cutoff_value * farsca + farloc

        return(list(farness_probs = probs, cutoff_value = cutoff_value[[1]]))
    } else {
        return(farness_probs = probs)
    }
}