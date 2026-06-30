#' Mallows Distance
#' 
#' Calculate the squared Mallows distance between all rows in data and the barycenter.
#' 
#' @details
#' The squared Mallows distance is defined according to the \code{LatentCase}:
#' \itemize{
#' \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{d_\mathrm{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#' \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{d_\mathrm{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}(\boldsymbol{r}-\boldsymbol{\mu}_R)
#'              +\mathbb{E}(U)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'  where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#' \item \code{"General"}: The latent variables do not have any nice properties:
#'    \deqn{d_\mathrm{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Delta}(\boldsymbol{r}-\boldsymbol{\mu}_R)
#'              +(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'  where:
#'  \itemize{
#'      \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'      \item \eqn{\boldsymbol{\Delta}=\text{diag}(\mathbb{E}(U^2_1),\dots,\mathbb{E}(U^2_p))/4}.
#'  }
#' }
#'
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @param mean_c (Optional) A vector specifying the mean of centers. Defaults to `NULL`, in which case 
#' it will be computed using the sample mean of centers.
#' @param mean_r (Optional) A vector specifying the mean of ranges Defaults to `NULL`, in which case 
#' it will be computed using the sample mean of ranges.
#' @return A vector with the squared Mallows distance of each observation.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_dist <- Mallows_dist(credit_card_int)
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
Mallows_dist <- function(data,mean_c=NULL,mean_r=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase

    if (is.null(mean_c) || is.null(mean_r)) {
        mean_c <- colMeans(C)
        mean_r <- colMeans(R)
    }

    C_0 <- t(C)-mean_c
    R_0 <- t(R)-mean_r
    
    if (case=="U_id_symmetric"){
        delta <- param[[1]]
        d2 <- diag(crossprod(C_0)) + delta*diag(crossprod(R_0))
    }else if (case=="U_id"){
        delta <- param[[1]]
        U_mean <- param[[2]]
        d2 <- diag(crossprod(C_0)) + 
                + delta*diag(crossprod(R_0)) +
                + U_mean*diag(crossprod(C_0,R_0))
    }else if (case=="General"){
        delta <- 1/4*diag(diag(param[[1]]))
        psi <- param[[2]]
        d2 <- diag(crossprod(C_0)) + 
                + diag(crossprod(R_0,delta%*%R_0)) +
                + diag(crossprod(C_0,psi%*%R_0))
    }
    names(d2) <- rownames(data)
    return(d2)
}

#' Interval-Mahalanobis Distance
#' 
#' Calculate the squared Interval-Mahalanobis distance of all rows in the data and the barycenter.
#' 
#' @details
#' The squared Interval-Mahalanobis distance between \eqn{\boldsymbol{x}=(\boldsymbol{c}^\top,\boldsymbol{r}^\top)^\top} and the barycenter \eqn{\boldsymbol{\mu}_B=(\boldsymbol{\mu}_C^\top,\boldsymbol{\mu}_R^\top)^\top} of a population with symbolic covariance matrix \eqn{\boldsymbol{\Sigma}_B} (see \code{\link{int_cov}}) is defined according to the \code{LatentCase}:
#' \itemize{
#' \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{d_\mathrm{IMah}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#' \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{\begin{aligned}
#'              d_\mathrm{IMah}(\boldsymbol{x})^2&=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\\
#'              &\quad+\mathbb{E}(U)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R),
#'     \end{aligned}}
#'  where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#' \item \code{"General"}: The latent variables do not have any nice properties:
#'    \deqn{\begin{aligned}
#'              d_\mathrm{IMah}(\boldsymbol{x})^2&=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\dfrac{1}{4}(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_{B}^{-1}\right)(\boldsymbol{r}-\boldsymbol{\mu}_R)\\
#'              &\quad+(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R),
#'    \end{aligned}}
#'  where:
#'  \itemize{
#'      \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt},
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'      \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'  }
#' }
#'
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @param z (Optional) A vector of 0 and 1, indicating which observations should be considered for the calculation.
#' If `z` is not `NULL`, `mean_c`, `mean_r`, and `cov` will be computed using only the observations with `z=1` (see \code{\link{int_mean_z}} and \code{\link{int_cov_z}}). 
#' Defaults to `NULL`.
#' @param mean_c (Optional) A vector specifying the mean of centers. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function, if `z` is also `NULL`.
#' @param mean_r (Optional) A vector specifying the mean of ranges. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function, if `z` is also `NULL`.
#' @param cov (Optional) A covariance matrix. Defaults to `NULL`, in which case it will be computed 
#' using the \code{\link{IMCD}} function, if `z` is also `NULL`.
#' @return A vector with the squared Interval-Mahalanobis distance of each observation.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute squared Interval-Mahalanobis distance using IMCD estimates of mean and covariance
#' credit_card_dist <- IMah_dist(credit_card_int)
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
IMah_dist <- function(data,z=NULL,mean_c=NULL,mean_r=NULL,cov=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase

    if(!is.null(z)) {
        mean_c <- int_mean_z(z, C)
        mean_r <- int_mean_z(z, R)
        cov <- int_cov_z(z,data)
    }
    
    if(is.null(mean_c) || is.null(mean_r) || is.null(cov)){
        IMCD_res <- IMCD(data,m=floor(data@NObs*0.75))
        if(is.null(cov)) cov <- IMCD_res$cov_IMCD
        if(is.null(mean_c)) mean_c <- IMCD_res$mean_IMCD_c
        if(is.null(mean_r)) mean_r <- IMCD_res$mean_IMCD_r
    }

    cov_inv <- safe_solve_cov(cov)
    
    if (case=="U_id_symmetric"){
        delta <- param[[1]]
        d2 <- mahalanobis(C, mean_c, cov_inv, inverted = TRUE) + 
                + delta*mahalanobis(R, mean_r, cov_inv, inverted = TRUE)
    }else if (case=="U_id"){
        delta <- param[[1]]
        U_mean <- param[[2]]
        C_0 <- t(C)-mean_c
        R_0 <- t(R)-mean_r
        d2 <- mahalanobis(C, mean_c, cov_inv, inverted = TRUE) + 
                + delta*mahalanobis(R, mean_r, cov_inv, inverted = TRUE) +
                + U_mean*diag(crossprod(C_0,cov_inv%*%R_0))
    }else if (case=="General"){
        e_UU <- param[[1]]
        psi <- param[[2]]
        C_0 <- t(C)-mean_c
        R_0 <- t(R)-mean_r
        d2 <- mahalanobis(C, mean_c, cov) + 
                + 1/4*diag(crossprod(R_0,(e_UU*cov_inv)%*%R_0)) +
                + diag(crossprod(C_0,cov_inv%*%psi%*%R_0))
    }
    names(d2) <- rownames(data)
    return(d2)
}

#' Interval-Mahalanobis distance for all pairs
#'
#' Calculate the squared Interval-Mahalanobis distance of all pairs of observations in the data.
#' 
#' @details
#' The squared Interval-Mahalanobis distance between \eqn{\boldsymbol{x}_1=(\boldsymbol{c}_1^\top,\boldsymbol{r}_1^\top)^\top} and \eqn{\boldsymbol{x}_2=(\boldsymbol{c}_2^\top,\boldsymbol{r}_2^\top)^\top} of a population with symbolic covariance matrix \eqn{\boldsymbol{\Sigma}_B} (see \code{\link{int_cov}}) is defined according to the \code{LatentCase}:
#' \itemize{
#'      \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'          \deqn{d_\mathrm{IMah}(\boldsymbol{x}_1,\boldsymbol{x}_2)^2=(\boldsymbol{c}_1-\boldsymbol{c}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_1-\boldsymbol{c}_2)+\delta(\boldsymbol{r}_1-\boldsymbol{r}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_1-\boldsymbol{r}_2),}
#'          where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'      \item \code{"U_id"}: The latent variables are identically distributed:
#'          \deqn{\begin{aligned}
#'              d_\mathrm{IMah}(\boldsymbol{x}_1,\boldsymbol{x}_2)^2&=(\boldsymbol{c}_1-\boldsymbol{c}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_1-\boldsymbol{c}_2)+\delta(\boldsymbol{r}_1-\boldsymbol{r}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_1-\boldsymbol{r}_2)\\
#'              &\quad+\mathbb{E}(U)(\boldsymbol{c}_1-\boldsymbol{c}_2)^\top\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_1-\boldsymbol{r}_2),
#'          \end{aligned}}
#'          where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#'      \item \code{"General"}: The latent variables do not have any nice properties:
#'          \deqn{\begin{aligned}
#'              d_\mathrm{IMah}(\boldsymbol{x}_1,\boldsymbol{x}_2)^2&=(\boldsymbol{c}_1-\boldsymbol{c}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_1-\boldsymbol{c}_2)+\dfrac{1}{4}(\boldsymbol{r}_1-\boldsymbol{r}_2)^{\top}\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_{B}^{-1}\right)(\boldsymbol{r}_1-\boldsymbol{r}_2)\\
#'              &\quad+(\boldsymbol{c}_1-\boldsymbol{c}_2)^{\top}\boldsymbol{\Sigma}_{B}^{-1}\boldsymbol{\Psi}(\boldsymbol{r}_1-\boldsymbol{r}_2),
#'          \end{aligned}}
#'          where:
#'          \itemize{
#'              \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'              \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'          }
#' }
#' If \code{cov} is not provided, it will be computed using the \code{\link{IMCD}} function.
#' Additionally, if \code{cov} is set as the identity matrix, the computed distance is the Mallows distance between pairs of observations.    
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @param cov (Optional) A covariance matrix. Defaults to `NULL`, in which case it will be computed 
#' using the \code{\link{IMCD}} function.
#' 
#' @return A matrix with the squared Interval-Mahalanobis distance of each pair of observations.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_dist <- IMah_dist_pairs(credit_card_int)
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
IMah_dist_pairs <- function(data,cov=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase
    n <- data@NObs
    
    if(is.null(cov)) {
        cov <- IMCD(data)$cov_IMCD
    }
    cov_inv <- safe_solve_cov(cov)

    D <- matrix(0, nrow = n, ncol = n)

    for(i in 1:n){
        for(j in 1:n){
            c_diff <- C[i, ] - C[j, ]
            r_diff <- R[i, ] - R[j, ]

            if (case=="U_id_symmetric"){
                delta <- param[[1]]
                d2 <- crossprod(c_diff, cov_inv %*% c_diff) + 
                        + delta * crossprod(r_diff, cov_inv %*% r_diff)
            }else if (case=="U_id"){
                delta <- param[[1]]
                U_mean <- param[[2]]
                d2 <- crossprod(c_diff, cov_inv %*% c_diff) + 
                        + delta * crossprod(r_diff, cov_inv %*% r_diff) +
                        + U_mean * crossprod(c_diff, cov_inv %*% r_diff)
            }else if (case=="General"){
                e_UU <- param[[1]]
                psi <- param[[2]]
                d2 <- crossprod(c_diff, cov_inv %*% c_diff) + 
                        + 1/4*crossprod(r_diff,(e_UU*cov_inv)%*%r_diff) +
                        + crossprod(c_diff,cov_inv%*%psi%*%r_diff)
            }
            D[i,j] <- d2
        }
    }
    colnames(D) <- rownames(D) <- rownames(data)
    return(D)
}
