#' Perform single iteration of C-step
#' @keywords internal
#' @param z A vector of 0 and 1, indicating which observations should be considered for the calculation
#' @param m An integer specifying number of observations to use
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return A list of z, covariance, barycenter and robust distances
c_step <- function(z,m,data){

    d2 <- IMah_dist(data=data,z=z)
    updated_z <- rep(0,length(z))
    updated_z[order(d2)[1:m]] <- 1

    S <- int_cov_z(updated_z,data)
    mean_c <- int_mean_z(updated_z,as.matrix(data@Centers))
    mean_r <- int_mean_z(updated_z,as.matrix(data@Ranges))

    return(list(updated_z=updated_z, S=S, mean_c=mean_c, 
                mean_r=mean_r, robust_dist=d2))
}

#' Randomly draw a subset of observations
#' @keywords internal
#' @param m An integer specifying the number of observations to use
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return A vector representing an m-length subset of X
draw_z <- function(m,data){
    n <- data@NObs; p <- data@NIVar
    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    z <- rep(0, n)
    j <- sample(n, p + 1)
    z[j] <- 1
    S <- int_cov_z(z,data)
    while(det(S) == 0 && length(j) < n - 1){
        j <- sample((1:n)[-j], 1)
        z[j] <- 1
        S <- int_cov_z(z,data)      
    }
    mean_c <- int_mean_z(z, C)
    mean_r <- int_mean_z(z, R)
    d2 <- IMah_dist(data = data, mean_c = mean_c, mean_r = mean_r, cov = S)
    z <- rep(0, n)
    z[order(d2)[1:m]] <- 1
    return(z)
}

#' Iterate through C-step
#' @keywords internal
#' @param z A vector of 0 and 1, indicating which observations should be considered for the calculation
#' @param m An integer specifying number of observations to use
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @param it An optional integer specifying the number of C-steps to perform.
#' With `it = 0`, C-step will be performed until convergence
#' @return A list of z, covariance, barycenter and robust distances
step_it <- function(z, m, data, it = 0){
    res <- c_step(z, m, data)
    det_new <- det(res$S)
    z <- res$updated_z

    if(it==0){
        det_old <- Inf
        while(det_old > det_new && det_new > 0){
            res <- c_step(z, m, data)
            det_old <- det_new
            det_new <- det(res$S)
            z <- res$updated_z
        }
    } else{
        for (i in seq_len(it-1)){
            res <- c_step(z, m, data)
            z <- res$updated_z
        }
    }
    return(res)
}

#' Choose the 10 best estimates after iterating twice through initial sets
#' @keywords internal
#' @param z_all A 2D matrix where each row specifies a subset of observations
#' @param m An integer specifying number of observations to use
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return A list of z, covariance, barycenter and robust distances
pick10 <- function(z_all, m, data){
    res <- apply(z_all, 2, function(z) {
        step_it(z, m, data, it = 2)
    })
    z <- sapply(res, function(x) {
        x$updated_z
        })
    det_S <- sapply(res, function(x) {
        det(x$S)
        })
    S <- lapply(res, function(x) {
        x$S
        })
    mean_c <- sapply(res, function(x) {
        x$mean_c
        })
    mean_r <- sapply(res, function(x) {
        x$mean_r
    })
    return(list(z = z[, head(order(det_S), 10)], S = S[head(order(det_S), 10)], 
                mean_c = mean_c[, head(order(det_S), 10)], 
                mean_r = mean_r[, head(order(det_S), 10)]))
}

#' Obtain unweighted estimates for data with <= 600 observations
#' @keywords internal
#' @param m An integer specifying the number of observations to use
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return A list of estimated barycenter and symbolic covariance matrix
smallIMCD <- function(m, data){
    
    # Sample initial subsets
    z_all <- sapply(1:500, function(i) {
        draw_z(m, data)
    })
    
    # Iterate for 10 best subsets
    z10 <- pick10(z_all, m, data)$z
    res10 <- apply(z10, 2, function(z) {
        step_it(z, m, data)
    })
    S_det <- sapply(res10, function(x) {
        det(x$S)
        })
    res <- res10[[which.min(S_det)]]
    return(res)
}

#' Obtain unweighted estimates for data with > 600 observations
#' @keywords internal
#' @param m An integer specifying number of observations to use
#' @param p An integer specifying the number of columns in X
#' @param n An integer specifying the number of total observations
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return A list of estimated location and scatter
bigIMCD <- function(m, p, n, data){
    k <- min(5, ceiling(n / 300))
    n_merge  <- min(1500, n)
    i_start <- rep(n_merge %% k, k)
    n_sub <- floor(n_merge / k)
    m_sub <- floor(n_sub * m / n)
    m_merge <- floor(n_merge * m / n)

    # Set indices to partition data
    if(!(n_merge %% k)){
        i_start <- n_sub * (0:(k - 1)) + 1
        i_end <- i_start - 1 + n_merge / k
    } else{
        i_start[seq(n_merge %% k)] <- seq(n_merge %% k) - 1
        i_start <- n_sub * (0:(k - 1)) + i_start + 1
        i_end <- c(i_start[2:k] - 1, n_merge)            
    }
    
    # Permute data and obtain sub-partitions
    samp <- sample(n)
    ind_z <- lapply(1:k, function(i) {
        samp[i_start[i]:i_end[i]]
        })
    data_sub <- lapply(1:k, function(i) {
        data[ind_z[[i]], ]
    })

    # Pick 10 best in each sub-partition
    data_merge <- do.call('rbind', data_sub)
    z_sub <- matrix(0, nrow = n_merge, ncol = 10*k)
    row_start <- 1
    for (i in seq_len(k)){
        z_all <- sapply(1:(500 / k), function(ind) { # ind is not used
            draw_z(m_sub, data_sub[[i]])
        })
        res10 <- pick10(z_all, m_sub, data_sub[[i]])

        # Actual number of rows and columns returned
        n_rows <- nrow(res10$z)
        n_cols <- ncol(res10$z)

        # Compute target indices based on actual block sizes
        row_idx <- row_start:(row_start + n_rows - 1)
        col_idx <- ((i - 1) * 10 + 1):((i - 1) * 10 + n_cols)

        # Safety: ensure we don’t exceed preallocated z_sub size
        row_idx <- row_idx[row_idx <= n_merge]
        col_idx <- col_idx[col_idx <= ncol(z_sub)]

        # Assign safely (truncate or pad if needed)
        if (length(row_idx) == nrow(res10$z)) {
            z_sub[row_idx, col_idx] <- res10$z
        } else {
            z_sub[row_idx, col_idx] <- res10$z[seq_along(row_idx), seq_along(col_idx), drop = FALSE]
        }

        row_start <- max(row_idx) + 1
    }
    # Perform 2 C-steps and select 10 best
    res10_z <- pick10(z_sub, m_merge, data_merge)$z
    z_merge <- apply(res10_z, 2, function(z) {
        as.numeric(data@ObsNames %in% data_merge[z==1,]@ObsNames)
    })
    
    # C-step until convergence for best subsets
    res <- NULL
    for (i in 1:10){
        res <- append(res, list(step_it(z_merge[, i], m, data)))
    }
    imin <- which.min(sapply(1:10, function(i) det(res[[i]]$S)))
    return(res[[imin]])
}

#' Interval-valued data Minimum Covariance Determinant (IMCD) estimation
#'
#' Applies an adaptation of the FAST-MCD algorithm to estimate location and scatter for interval-valued data.
#' 
#' @param data An \code{\linkS4class{intData}} object containing the interval-valued dataset (macrodata).
#' @param m An integer specifying the subset size to use for the estimation. Defaults to `floor(0.75*n)`.
#' @param cutoff Indicates which cutoff should be considered for reweighting the estimates:
#' \itemize{
#'    \item \code{"chi-squared"}: The traditional 97.5\% Chi-Squared quantile.
#'    \item \code{"raw"}: No reweighting.
#'    \item \code{"adjbox"}: Adjusted Boxplots (package \code{robustbase}).
#'    \item \code{"F-dist"}: The quantile of the scaled F distribution (adapted from package \code{CerioliOutlierDetection}).
#'    \item \code{"farness"}: "Farness" is estimated from the robust distance (adapted from package \code{cellWise}).
#' }
#' Defaults to \code{"farness"}.
#' @param cutoff_lvl A numeric value specifying the level of the cutoff to be used. 
#' \itemize{
#'      \item If \code{cutoff="chi-squared"}, \code{cutoff_lvl} is the quantile of the Chi-squared distribution (default is 0.975).
#'      \item If \code{cutoff="adjbox"}, \code{cutoff_lvl} is the coefficient for the adjusted boxplot (default is 1.5).
#'      \item If \code{cutoff="F-dist"}, \code{cutoff_lvl} is the quantile of the F-distribution (default is 0.975).
#'      \item If \code{cutoff="farness"}, \code{cutoff_lvl} represents the threshold for farness, with a default of 0.99.
#'      \item If \code{cutoff="raw"}, \code{cutoff_lvl} is ignored.
#' }
#' If no value is provided, the function uses the default values associated with each cutoff method.
#' @return A list containing the robustly estimated parameters:
#'   \item{\code{mean_IMCD_c}}{Estimated mean of the centers of the interval data.}
#'   \item{\code{mean_IMCD_r}}{Estimated mean of the ranges of the interval data.}
#'   \item{\code{cov_IMCD}}{Estimated covariance (scatter) matrix (\code{\link{int_cov}}) for the data.}
#'   \item{\code{final_z}}{Binary vector indicating the inclusion of each observation in the reweighted subset.}
#'   \item{\code{cutoff}}{The cutoff method used for reweighting.}
#'   \item{\code{cutoff_value}}{Cutoff value used for reweighting.}
#'   \item{\code{robust_dist}}{Robust distances (\code{\link{IMah_dist}}) for each observation.}
#'   \item{\code{farness_probs}}{Farness probabilities (if \code{cutoff} is set to \code{"farness"}).}
#' @export
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
#' @references Adapted from \url{https://github.com/frankp-0/fastMCD}.
#' @references The case \code{cutoff=="F-dist"} is adapted from package \code{CerioliOutlierDetection} (\url{https://cran.r-project.org/package=CerioliOutlierDetection}).
#' @examples
#' # Example using creditcard dataset
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Obtain reweighted IMCD estimates using farness cutoff
#' credit_card_IMCD <- IMCD(credit_card_int, 
#'                          m = floor(nrow(credit_card_int)*0.75), 
#'                          cutoff = "farness", 
#'                          cutoff_lvl = 0.9)
IMCD <- function(data, 
                m = 0, 
                cutoff=c("farness","adjbox","chi-squared","F-dist","raw"), 
                cutoff_lvl=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")
    
    cutoff <- match.arg(cutoff)
    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    n <- data@NObs; p <- data@NIVar

    if(p==1){stop("data needs to have at least 2 variables.")}

    if (is.null(cutoff_lvl)){
        cutoff_lvl <- switch(cutoff,
                            "chi-squared" = 0.975,
                            "adjbox" = 1.5,
                            "F-dist" = 0.975,
                            "farness" = 0.95,
                            0.975)
    }
    
    # Set m
    if (!m) {
        m <- floor(0.75*n)
    }
    
    # Call helper for big or small data
    if (m == n) {
        z <- rep(1, n)
        d2 <- IMah_dist(data=data,z=z)
    } else {
        if (n <= 600) {
            res <- smallIMCD(m, data)
        } else {
            res <- bigIMCD(m, p, n, data)
        }                        
        z <- res$updated_z
        d2 <- res$robust_dist
    }
    
    farness_probs <- NA

    # Reweight estimates
    if (cutoff=="chi-squared"){
        cutoff_value <- qchisq(cutoff_lvl, df = p)
        w <- ifelse(d2 <= cutoff_value, 1, 0)
    }else if (cutoff=="raw"){
        cutoff_value <- NA
        w <- z
    }else if (cutoff=="adjbox"){
        if (!requireNamespace("robustbase", quietly = TRUE)) {
            stop("Package 'robustbase' is required for cutoff=='adjbox'.")
        }
        cutoff_value <- robustbase::adjboxStats(d2, coef=cutoff_lvl, doScale = FALSE)$fence
        w <- ifelse((d2 >= cutoff_value[1])&(d2 <= cutoff_value[2]), 1, 0)
    }else if (cutoff=="F-dist"){
        if (!requireNamespace("CerioliOutlierDetection", quietly = TRUE)) {
            stop("Package 'CerioliOutlierDetection' is required for cutoff=='F-dist'.")
        }
        delta <- 1-cutoff_lvl
        hr05 <- CerioliOutlierDetection::hr05CutoffMvnormal(n, p, m/n, delta)
        dfz <- hr05$m.pred - p + 1
        cutoff_value <- hr05$m.pred * p * qf(1 - delta, df1 = p, df2 = dfz) / dfz
        w <- ifelse(d2 <= cutoff_value, 1, 0)
    }else if(cutoff=="farness"){
        farness_results <- farness(d2, cutoff_value = cutoff_lvl)
        farness_probs <- farness_results$farness_probs
        cutoff_value <- farness_results$cutoff_value
        w <- ifelse(farness_probs <= cutoff_lvl, 1, 0)
    }

    reweighted_mu_c <- int_mean_z(w, C)
    reweighted_mu_r <- int_mean_z(w, R)
    reweighted_sigma <- int_cov_z(w, data)
    reweighted_d2 <- IMah_dist(data,mean_c=reweighted_mu_c,mean_r=reweighted_mu_r,cov=reweighted_sigma)
    names(w) <- rownames(data)

    return(list("mean_IMCD_c"=reweighted_mu_c, "mean_IMCD_r"=reweighted_mu_r, 
                "cov_IMCD"=reweighted_sigma, "final_z"=w, "cutoff"=cutoff,
                "cutoff_value"=cutoff_value, "robust_dist"=reweighted_d2, "farness_probs"=farness_probs))
}

