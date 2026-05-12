#' Relative Frobenius Error
#' 
#' Computes the relative Frobenius error between an estimated covariance matrix and the ground truth.
#' 
#' @details 
#' The relative Frobenius error is given by: 
#' \deqn{\dfrac{\|\boldsymbol{A} - \boldsymbol{B}\|_F}{\|\boldsymbol{B}\|_F}=\dfrac{\sqrt{\sum\limits_{i=1}^{p}\sum\limits_{j=1}^{p}|[\boldsymbol{A}]_{ij}-[\boldsymbol{B}]_{ij}|^2}}{\sqrt{\sum\limits_{i=1}^{p}\sum\limits_{j=1}^{p}|[\boldsymbol{B}]_{ij}|^2}},}
#' where \eqn{\boldsymbol{A}} and \eqn{\boldsymbol{B}} are the estimated and ground truth covariance matrices, respectively.
#' 
#' @param est_cov Estimated covariance matrix.
#' @param ground_truth_cov Ground truth covariance matrix.
#' @return Frobenius error between the two matrices.
#' @export
frobenius_error <- function(est_cov, ground_truth_cov) {
  return(sqrt(sum((est_cov - ground_truth_cov)^2))/sqrt(sum((ground_truth_cov)^2)))
}

#' Angle Error
#' 
#' Computes the angle error between eigenvalues of the estimated covariance matrix and of the ground truth covariance matrix.
#' 
#' @details
#' The angle error is given by:
#' \deqn{1-\dfrac{\hat{\boldsymbol{a}}^\top\boldsymbol{a}}{\sqrt{\hat{\boldsymbol{a}}^\top\hat{\boldsymbol{a}}}\sqrt{\boldsymbol{a}^\top\boldsymbol{a}}},}
#' where \eqn{\hat{\boldsymbol{a}}} and \eqn{\boldsymbol{a}} are the eigenvalues of the estimated and ground truth covariance matrices, respectively.
#' 
#' @param est_cov Estimated covariance matrix.
#' @param ground_truth_cov Ground truth covariance matrix.
#' @return Angle error between eigenvalues.
#' @importFrom geigen geigen
#' @export
angle_error <- function(est_cov, ground_truth_cov){
    a_est <- sort(geigen::geigen(est_cov, diag(nrow(est_cov)), only.values=TRUE)$values)
    a <- sort(geigen::geigen(ground_truth_cov, diag(nrow(ground_truth_cov)), only.values=TRUE)$values)

    return(as.numeric(1-(t(a_est)%*%a)/(sqrt(t(a_est)%*%a_est)*sqrt(t(a)%*%a))))
}

#' Kullback-Leibler (KL) Divergence
#' 
#' Computes the Kullback-Leibler (KL) divergence between an estimated covariance matrix and the ground truth. Assumes normal multivariate distributions.
#' 
#' @details
#' The KL divergence between two \eqn{p}-dimensional Gaussians \eqn{\mathcal{N}(\boldsymbol{\mu}, \hat{\boldsymbol{\Sigma}})} and \eqn{\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})} is given by:
#' \deqn{\dfrac{1}{2}\left(\text{tr}(\boldsymbol{\Sigma}^{-1}\hat{\boldsymbol{\Sigma}}) + \log\left(\dfrac{\det(\boldsymbol{\Sigma})}{\det(\hat{\boldsymbol{\Sigma}})}\right) - p\right),}
#' where \eqn{\hat{\boldsymbol{\Sigma}}} and \eqn{\boldsymbol{\Sigma}} are the estimated and ground truth covariance matrices, respectively.
#' 
#' @param est_cov Estimated covariance matrix.
#' @param ground_truth_cov Ground truth covariance matrix.
#' @return KL divergence between the two matrices.
#' @references Yufeng Zhang, Wanwei Liu, Zhenbang Chen, Ji Wang, and Kenli Li. On the properties of Kullback-Leibler divergence between multivariate gaussian distributions, 2023. \url{https://arxiv.org/abs/2102.05485}
#' @export
KL_divergence <- function(est_cov, ground_truth_cov) {
    if (det(est_cov) <= 0 || det(ground_truth_cov) <= 0) {
        stop("Covariance matrices must be positive definite.")
    }
    S <- est_cov %*% solve(ground_truth_cov)
    return(sum(diag(S))-log(det(S))-nrow(ground_truth_cov))
}