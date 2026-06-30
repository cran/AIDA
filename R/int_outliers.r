#' Outlier Detection for Interval-Valued Data Based on Robust Distances
#'
#' Identifies potential outliers in interval-valued data using robust distance-based methods with customizable cutoff criteria.
#'
#' This function classifies observations as outliers based on robust distances and user-defined cutoff methods. It supports various approaches, including Chi-Squared quantiles, adjusted boxplots, F distribution quantiles, and farness probabilities.
#'
#' @param robust_dist A numeric vector containing the robust distances for each observation.
#' @param cutoff_lvl A numeric value specifying the level of the cutoff to be used. 
#' \itemize{
#'      \item If \code{cutoff="chi-squared"}, \code{cutoff_lvl} is the quantile of the Chi-squared distribution (default is 0.975).
#'      \item If \code{cutoff="adjbox"}, \code{cutoff_lvl} is the coefficient for the adjusted boxplot (default is 1.5).
#'      \item If \code{cutoff="F-dist"}, \code{cutoff_lvl} is the quantile of the F-distribution (default is 0.975).
#'      \item If \code{cutoff="farness"}, \code{cutoff_lvl} represents the threshold for farness, with a default of 0.99.
#'      \item If \code{cutoff="raw"}, \code{cutoff_lvl} is ignored.
#' }
#' If no value is provided, the function uses the default values associated with each cutoff method.
#' @param cutoff A character string specifying the method for setting the outlier cutoff threshold. Options include:
#'   \itemize{
#'      \item \code{"chi-squared"}: Outliers are identified based on a specified Chi-Squared quantile.
#'      \item \code{"adjbox"}: Uses adjusted boxplot statistics (from \code{robustbase}) to classify outliers.
#'      \item \code{"F-dist"}: Applies a cutoff derived from the F and Beta distributions for robust outlier detection.
#'      \item \code{"farness"}: Identifies outliers based on a "farness" threshold, determined by the robust distance distribution.
#'   }
#'   Default is \code{"farness"}.
#' @param cutoff_lvl A numeric value specifying the level of the cutoff to be used. 
#' \itemize{
#'      \item If \code{cutoff="chi-squared"}, \code{cutoff_lvl} is the quantile of the Chi-squared distribution (default is 0.975).
#'      \item If \code{cutoff="adjbox"}, \code{cutoff_lvl} is the coefficient for the adjusted boxplot (default is 1.5).
#'      \item If \code{cutoff="F-dist"}, \code{cutoff_lvl} is the significance level for identifying outliers (default is 0.95).
#'      \item If \code{cutoff="farness"}, \code{cutoff_lvl} represents the threshold for farness, with a default of 0.99.
#' }
#' If no value is provided, the function uses the default values associated with each cutoff method.
#' @param p The number of variables in the data. Required for \code{"chi-squared"} and \code{"F-dist"} cutoff methods.
#' @param z A binary vector indicating the subset of observations used for initial robust estimation. Required for the \code{"F-dist"} cutoff method.
#' 
#' @return A list with the following components:
#'   \item{\code{outliers_names}}{Character vector of names for observations classified as outliers.}
#'   \item{\code{is_outlier}}{Logical vector indicating whether each observation is an outlier (TRUE) or not (FALSE).}
#'   \item{\code{cutoff}}{The cutoff method used for detecting outliers.}
#'   \item{\code{cutoff_value}}{Cutoff value used for detecting outliers.}
#'   \item{\code{farness_probs}}{Numeric vector of farness probabilities for each observation (only if \code{cutoff} is set to \code{"farness"}).}
#' 
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
#' @references Case \code{cutoff=="F-dist"} is adapted from package \code{CerioliOutlierDetection} (\url{https://cran.r-project.org/package=CerioliOutlierDetection}).
#' @examples
#' # Example of detecting outliers using robust distances
#' set.seed(42)
#' robust_dist <- abs(rnorm(100))
#' result <- int_outliers(robust_dist, cutoff = "chi-squared", p = 5)
#' 
#' # Example using creditcard dataset
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute robust distances using IMCD estimates of mean and covariance
#' credit_card_dist <- IMah_dist(credit_card_int)
#' 
#' # Detect outliers using farness cutoff
#' credit_card_outliers <- int_outliers(credit_card_dist, 
#'                                      cutoff = "farness", 
#'                                      cutoff_lvl = 0.9)
#' @export
int_outliers <- function(robust_dist,
                        cutoff=c("farness","adjbox","chi-squared","F-dist"),
                        cutoff_lvl=NULL,
                        p=NULL,
                        z=NULL){
    cutoff <- match.arg(cutoff)

    if (is.null(cutoff_lvl)){
        cutoff_lvl <- switch(cutoff,
                            "chi-squared" = 0.975,
                            "adjbox" = 1.5,
                            "F-dist" = 0.95,
                            "farness" = 0.95,
                            0.975)
    }

    if((cutoff=="chi-squared"||cutoff=="F-dist")&&is.null(p)){stop("For cutoff='chi-squared' or cutoff='F-dist', you must provide the number of variables p.")}

    if(cutoff=="F-dist"&&is.null(z)){stop("For cutoff='F-dist', you must provide z.")}

    farness_probs <- NA
    if (cutoff=="chi-squared"){
        cutoff_value <- qchisq(cutoff_lvl, df = p)
        w <- ifelse(robust_dist <= cutoff_value, FALSE, TRUE)
    }else if (cutoff=="adjbox"){
        if (!requireNamespace("robustbase", quietly = TRUE)) {
            stop("Package 'robustbase' is required for cutoff=='adjbox'.")
        }
        cutoff_value <- robustbase::adjboxStats(robust_dist, coef=cutoff_lvl, doScale = FALSE)$fence
        w <- ifelse((robust_dist >= cutoff_value[1])&(robust_dist <= cutoff_value[2]), FALSE, TRUE)
    }else if (cutoff=="F-dist"){
        critfcn <- function(mm, vv, ww) {
            function(siga) {
                if (mm < vv) 
                stop("DF2 parameter for F distribution is negative.")
                if (mm <= vv + 1)
                stop("Shape2 parameter for beta distribution will be zero or negative")
                if (any(is.na(ww)))
                stop("There are missing weights.")

                crit.w0 <- ( (mm * mm - 1) * vv ) * qf(1 - siga, df1 = vv, df2 = mm - vv) / (mm * (mm - vv))
                crit.w1 <- (mm - 1) * (mm - 1) * qbeta(1 - siga, shape1 = vv / 2, shape2 = (mm - vv - 1) / 2) / mm

                critval <- rep(NA, length(ww))
                critval[ww == 1] <- crit.w1 
                critval[ww == 0] <- crit.w0
                critval
            }
        }
        critvalfcn <- critfcn(sum(z), p, z)
        signif.alpha <- 1-cutoff_lvl
        cutoff_value <- sapply(signif.alpha, critvalfcn)
        w <- (robust_dist > cutoff_value)[,1]
    }else if(cutoff=="farness"){
        farness_results <- farness(robust_dist, cutoff_value = cutoff_lvl)
        farness_probs <- farness_results$farness_probs
        cutoff_value <- farness_results$cutoff_value
        w <- ifelse(farness_probs <= cutoff_lvl, FALSE, TRUE)
        
    }
    names(w) <- names(robust_dist)
    outliers <- names(which(w))
    return(list("outliers_names"=outliers,"is_outlier"=w,"cutoff"=cutoff,"cutoff_value"=cutoff_value,"farness_probs"=farness_probs))
}