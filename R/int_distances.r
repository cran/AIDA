#' Mallows Distance
#' 
#' Calculate the squared Mallows distance between all rows in data and the barycenter.
#' 
#' @details
#' The squared Mallows distance is defined according to the \code{LatentCase}:
#' \itemize{
#' \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{d_{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#' \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{d_{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}(\boldsymbol{r}-\boldsymbol{\mu}_R)
#'              +\mathbb{E}(U)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'  where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#' \item \code{"General"}: The latent variables do not have any nice properties:
#'    \deqn{d_{M}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}(\boldsymbol{c}-\boldsymbol{\mu}_C)+(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Delta}(\boldsymbol{r}-\boldsymbol{\mu}_R)
#'              +(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'  where:
#'  \itemize{
#'      \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'      \item \eqn{\boldsymbol{\Delta}=\text{diag}(\mathbb{E}(U^2_1),\dots,\mathbb{E}(U^2_p))/4}.
#'  }
#' }
#'
#' @param data An \linkS4class{intData} object containing the macrodata/interval data
#' @param mean_c The mean vector of the centers
#' @param mean_r The mean vector of the ranges
#' @return A vector with the squared Mallows distance of each observation.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_dist<-Mallows_dist(credit_card_int)
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
Mallows_dist <- function(data,mean_c=NULL,mean_r=NULL){

    C<-as.matrix(data@Centers)
    R<-as.matrix(data@Ranges)
    param<-data@LatentParam
    case<-data@LatentCase

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
    names(d2)<-rownames(data)
    return(d2)
}

#' Interval-Mahalanobis Distance
#' 
#' Calculate the squared Interval-Mahalanobis distance of all rows in the data and the barycenter.
#' 
#' @details
#' The squared Interval-Mahalanobis distance is defined according to the \code{LatentCase}:
#' \itemize{
#' \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{d_{IMah}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R),}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#' \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{d_{IMah}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\\
#'              +\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)+\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C),}
#'  where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#' \item \code{"General"}: The latent variables do not have any nice properties:
#'    \deqn{d_{IMah}(\boldsymbol{x})^2=(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)+\dfrac{1}{4}(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_{B}^{-1}\right)(\boldsymbol{r}-\boldsymbol{\mu}_R)\\
#'              +\dfrac{1}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)^{\top}\boldsymbol{\Sigma}_{B}^{-1}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R)+\dfrac{1}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)^{\top}\boldsymbol{\Psi}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C),}
#'  where:
#'  \itemize{
#'      \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)}, \eqn{i\neq j}, with \eqn{\mathcal{E}(U_i,U_j)=\int_0^1 F_{U_i}^{-1}(t) F_{U_j}^{-1}(t) \, dt},
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ii}=\mathbb{E}(U_i^2)}, \eqn{i,j=1,\dots,p},
#'      \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'  }
#' }
#'
#' @param data An \linkS4class{intData} object containing the macrodata/interval data
#' @param z A vector of 0 and 1, indicating which observations should be considered for the calculation. 
#' You must provide either \code{z} or (\code{mean_c}, \code{mean_r} and \code{cov})
#' @param mean_c The mean vector of the centers
#' @param mean_r The mean vector of the ranges
#' @param cov The symbolic covariance matrix
#' @return A vector with the squared Interval-Mahalanobis distance of each observation.
#' @importFrom stats mahalanobis
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' z <- rep(1, nrow(credit_card_int))
#' credit_card_dist<-IMah_dist(credit_card_int,z)
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
IMah_dist <- function(data,z=NULL,mean_c=NULL,mean_r=NULL,cov=NULL){

    if (is.null(z) && (is.null(mean_c) || is.null(mean_r) || is.null(cov))) {
        stop("You must provide either z or (mean_c, mean_r and cov)")
    }

    C<-as.matrix(data@Centers)
    R<-as.matrix(data@Ranges)
    param<-data@LatentParam
    case<-data@LatentCase

    if (is.null(mean_c) || is.null(mean_r) || is.null(cov)) {
        mean_c <- int_mean_z(z, C)
        mean_r <- int_mean_z(z, R)
        cov <- int_cov_z(z,data)
    }

    cov_inv <- solve(cov)
    
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
                + U_mean/2*diag(crossprod(C_0,cov_inv%*%R_0)) + 
                + U_mean/2*diag(crossprod(R_0,cov_inv%*%C_0))
    }else if (case=="General"){
        e_UU <- param[[1]]
        psi <- param[[2]]
        C_0 <- t(C)-mean_c
        R_0 <- t(R)-mean_r
        d2 <- mahalanobis(C, mean_c, cov) + 
                + 1/4*diag(crossprod(R_0,(e_UU*cov_inv)%*%R_0)) +
                + 1/2*diag(crossprod(C_0,cov_inv%*%psi%*%R_0)) + 
                + 1/2*diag(crossprod(R_0,psi%*%cov_inv%*%C_0))
    }
    names(d2)<-rownames(data)
    return(d2)
}

#' Interval-Mahalanobis distance for all pairs
#'
#' Calculate the squared Interval-Mahalanobis distance of all pairs of observations in the data.
#' 
#' @details
#' The squared Interval-Mahalanobis distance is defined according to the \code{LatentCase}:
#' \itemize{
#'      \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'          \deqn{d_{IMah}(\boldsymbol{x}_i,\boldsymbol{x}_j)^2=(\boldsymbol{c}_i-\boldsymbol{c}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_i-\boldsymbol{c}_j)+\delta(\boldsymbol{r}_i-\boldsymbol{r}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_i-\boldsymbol{r}_j),}
#'          where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'      \item \code{"U_id"}: The latent variables are identically distributed:
#'          \deqn{d_{IMah}(\boldsymbol{x}_i,\boldsymbol{x}_j)^2=(\boldsymbol{c}_i-\boldsymbol{c}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_i-\boldsymbol{c}_j)+\delta(\boldsymbol{r}_i-\boldsymbol{r}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_i-\boldsymbol{r}_j)\\
#'              +\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{c}_i-\boldsymbol{c}_j)^\top\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{r}_i-\boldsymbol{r}_j)+\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{r}_i-\boldsymbol{r}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_i-\boldsymbol{c}_j),}
#'          where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#'      \item \code{"General"}: The latent variables do not have any nice properties:
#'          \deqn{d_{IMah}(\boldsymbol{x}_i,\boldsymbol{x}_j)^2=(\boldsymbol{c}_i-\boldsymbol{c}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_i-\boldsymbol{c}_j)+\dfrac{1}{4}(\boldsymbol{r}_i-\boldsymbol{r}_j)^{\top}\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_{B}^{-1}\right)(\boldsymbol{r}_i-\boldsymbol{r}_j)\\
#'              +\dfrac{1}{2}(\boldsymbol{c}_i-\boldsymbol{c}_j)^{\top}\boldsymbol{\Sigma}_{B}^{-1}\boldsymbol{\Psi}(\boldsymbol{r}_i-\boldsymbol{r}_j)+\dfrac{1}{2}(\boldsymbol{r}_i-\boldsymbol{r}_j)^{\top}\boldsymbol{\Psi}\boldsymbol{\Sigma}_{B}^{-1}(\boldsymbol{c}_i-\boldsymbol{c}_j),}
#'          where:
#'          \itemize{
#'              \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)}, \eqn{i\neq j}, with \eqn{\mathcal{E}(U_i,U_j)=\int_0^1 F_{U_i}^{-1}(t) F_{U_j}^{-1}(t) \, dt},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ii}=\mathbb{E}(U_i^2)}, \eqn{i,j=1,\dots,p},
#'              \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'          }
#' }
#' @param data An \linkS4class{intData} object containing the macrodata/interval data
#' @param cov The symbolic covariance matrix
#' 
#' @return A matrix with the squared Interval-Mahalanobis distance of each pair of observations.
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_dist<-IMah_dist_pairs(credit_card_int)
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
IMah_dist_pairs <- function(data,cov=NULL){

    C<-as.matrix(data@Centers)
    R<-as.matrix(data@Ranges)
    param<-data@LatentParam
    case<-data@LatentCase
    n <- data@NObs
    
    if(is.null(cov)) {
        cov <- IMCD(data)$cov_IMCD
    }
    cov_inv <- solve(cov)

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
                        + U_mean/2 * crossprod(c_diff, cov_inv %*% r_diff) + 
                        + U_mean/2 * crossprod(r_diff, cov_inv %*% c_diff)
            }else if (case=="General"){
                e_UU <- param[[1]]
                psi <- param[[2]]
                d2 <- crossprod(c_diff, cov_inv %*% c_diff) + 
                        + 1/4*crossprod(r_diff,(e_UU*cov_inv)%*%r_diff) +
                        + 1/2*crossprod(c_diff,cov_inv%*%psi%*%r_diff) + 
                        + 1/2*crossprod(r_diff,psi%*%cov_inv%*%c_diff)
            }
            D[i,j] <- d2
        }
    }
    names(D) <- rownames(D) <- rownames(data)
    return(D)
}
