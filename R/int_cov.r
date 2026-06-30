#' Sample Mean
#' 
#' Calculate the mean of X in function of z
#' 
#' @details 
#' This function calculates the mean of \eqn{\boldsymbol{X}} in function of \eqn{\boldsymbol{z}}. If \eqn{\boldsymbol{z}} is a vector of 0 and 1, the mean is calculated for the \eqn{m} observations that are equal to 1:
#' \deqn{\bar{\boldsymbol{x}}(\boldsymbol{z}) = \dfrac{1}{m} \boldsymbol{X}^\top \boldsymbol{z}.}
#' 
#' @param z A vector of 0 and 1, indicating which observations should be considered for the calculation
#' @param X A matrix where the rows correspond to observations and the columns to variables
#' @return A vector where each element is the mean for each variable
#' @export
#' @examples
#' n <- 100
#' p <- 4
#' X <- matrix(rnorm(n * p), ncol = p)
#' 
#' # if we consider all the observations the result obtained is the same as colMeans()
#' z <- c(rep(1, n))
#' int_mean_z(z, X)
#' colMeans(X)
int_mean_z <- function(z,X){
    return(c(1/sum(z)*crossprod(X,z)))
}

#' Interval-valued Covariance
#' 
#' Calculate the interval-valued covariance matrix based on the covariance matrices of the centers and ranges or data.
#' 
#' @details 
#' This function calculates the interval-valued covariance matrix, \eqn{\boldsymbol{\Sigma}_B}, based on the covariance matrices of the centers, \eqn{\boldsymbol{\Sigma}_{CC}}, ranges, \eqn{\boldsymbol{\Sigma}_{RR}}, and the covariance matrix between the centers and ranges, \eqn{\boldsymbol{\Sigma}_{CR}=\boldsymbol{\Sigma}_{RC}^\top}.
#' The covariance matrix is defined according to the \code{LatentCase}:
#' \itemize{
#'  \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{\boldsymbol{\Sigma}_B=\boldsymbol{\Sigma}_{CC}+\delta\boldsymbol{\Sigma}_{RR},}
#'    where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'  \item \code{"U_id"}: The latent variables are identically distributed:
#'      \deqn{\boldsymbol{\Sigma}_B=\boldsymbol{\Sigma}_{CC}+\delta\boldsymbol{\Sigma}_{RR}+\dfrac{\mathbb{E}(U)}{2}\left(\boldsymbol{\Sigma}_{CR}+\boldsymbol{\Sigma}_{RC}\right),}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameters of the latent variables.
#'  \item \code{"General"}: The latent variables do not have any nice properties:
#'      \deqn{\boldsymbol{\Sigma}_B=\boldsymbol{\Sigma}_{CC}+\dfrac{1}{4}\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_{RR}+\dfrac{1}{2}\boldsymbol{\Sigma}_{CR}\boldsymbol{\Psi}+\dfrac{1}{2}\boldsymbol{\Psi}\boldsymbol{\Sigma}_{RC}}
#'    where:
#'      \itemize{
#'          \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt}, 
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'          \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'     }
#' }
#' The covariance matrix can be calculated either based on the covariance matrices of the centers and ranges or based on the data. If the data is provided, the covariance matrices are calculated using the sample covariance of the centers and ranges and the sample covariance between centers and ranges.
#' For the robust estimation of the covariance matrix, see \code{\link{IMCD}}.
#' 
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data. 
#' If \code{data} is provided, the covariance matrix is calculated based on the the sample covariance of the centers and ranges and the sample covariance between centers and ranges, and the parameters of the latent variables contained in the \code{intData} object. 
#' If \code{data} is not provided, the covariance matrix is calculated based on `sigma_cc`, `sigma_rr`, `sigma_cr`, `LatentParam`, and `LatentCase`.
#' @param sigma_cc Covariance matrix of the centers.
#' @param sigma_rr Covariance matrix of the ranges.
#' @param sigma_cr Covariance matrix between the centers and ranges.
#' @param LatentParam A list with the parameters of the latent variables.
#'  Expects a list with a single number if `LatentCase` is `"U_id_symmetric"`, a list of two numbers if `LatentCase` is `"U_id"`, and a list of two matrices if `LatentCase` is `"General"`.
#' @param LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @return The symbolic covariance matrix.
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' credit_card_cov <- int_cov(credit_card_int)
int_cov <- function(data=NULL,
                    sigma_cc=NULL,
                    sigma_rr=NULL,
                    sigma_cr=NULL,
                    LatentParam=NULL,
                    LatentCase=c("U_id_symmetric","U_id","General")){
    
    if(!is.null(data)){
        if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")
        LatentParam <- data@LatentParam
        case <- data@LatentCase
        C <- as.matrix(data@Centers)
        R <- as.matrix(data@Ranges)
        sigma_cc <- cov(C)
        sigma_rr <- cov(R)
        sigma_cr <- cov(C,R)
    }else if(!is.null(sigma_cc)&&!is.null(sigma_rr)&&!is.null(LatentParam)){
        case <- match.arg(LatentCase)
        if(is.null(sigma_cr)&&case!="U_id_symmetric") stop("sigma_cr is missing.")
    }else{
       stop("Must provide either an intData object or the covariance matrices and the latent variables parameters and case.")
    }

    if (case=="U_id_symmetric"){
        delta <- LatentParam[[1]]
        sigma_b <- sigma_cc+delta*sigma_rr
    }else if (case=="U_id") {
        delta <- LatentParam[[1]]
        U_mean <- LatentParam[[2]]
        sigma_b <- sigma_cc+delta*sigma_rr+U_mean/2*(sigma_cr+t(sigma_cr))
    }else if (case=="General") {
        e_UU <- LatentParam[[1]]
        psi <- LatentParam[[2]]
        sigma_b <- sigma_cc+1/4*e_UU*sigma_rr+1/2*sigma_cr%*%psi+1/2*psi%*%t(sigma_cr)
    }
    if (!is.null(data)) {
        colnames(sigma_b) <- rownames(sigma_b) <- colnames(data)
    }
    return(sigma_b)
}

#' Sample Interval-valued Covariance
#' 
#' Calculate the interval-valued covariance matrix in function of z
#' 
#' @details 
#' Let \eqn{\boldsymbol{z}\in\{0,1\}^n} be a vector indicating which \eqn{m} observations are ``active''. This function calculates the sample interval-valued covariance matrix in function of \eqn{\boldsymbol{z}}: \eqn{\boldsymbol{S}_B(\boldsymbol{z})}.
#' Let \eqn{\boldsymbol{C}}, \eqn{\boldsymbol{R}} be the matrices of centers and ranges, respectively. Additionally, set:
#' \deqn{\overline{\boldsymbol{c}}_B(\boldsymbol{z})=\dfrac{1}{m}\boldsymbol{C}^{\top}\boldsymbol{z}, \qquad \overline{\boldsymbol{r}}_B(\boldsymbol{z})=\dfrac{1}{m}\boldsymbol{R}^{\top}\boldsymbol{z}.}
#' The sample interval-valued covariance matrix is obtained according to the \code{LatentCase}:
#' \itemize{
#'  \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{\boldsymbol{S}_B(\boldsymbol{z})=\left(\dfrac{1}{m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{c}_{h}\boldsymbol{c}_{h}^{\top}\right)-\overline{\boldsymbol{c}}_B(\boldsymbol{z})\overline{\boldsymbol{c}}_B(\boldsymbol{z})^\top+\left(\dfrac{\delta}{m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{r}_{h}\boldsymbol{r}_{h}^{\top}\right)-\delta\overline{\boldsymbol{r}}_B(\boldsymbol{z})\overline{\boldsymbol{r}}_B(\boldsymbol{z})^\top,}
#'    where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'  \item \code{"U_id"}: The latent variables are identically distributed:
#'      \deqn{\begin{aligned}
#'              \boldsymbol{S}_B(\boldsymbol{z})&=\left(\dfrac{1}{m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{c}_{h}\boldsymbol{c}_{h}^{\top}\right)-\overline{\boldsymbol{c}}_B(\boldsymbol{z})\overline{\boldsymbol{c}}_B(\boldsymbol{z})^\top+\left(\dfrac{\delta}{m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{r}_{h}\boldsymbol{r}_{h}^{\top}\right)-\delta\overline{\boldsymbol{r}}_B(\boldsymbol{z})\overline{\boldsymbol{r}}_B(\boldsymbol{z})^\top\\
#'              &\quad+\left(\dfrac{\mathbb{E}(U)}{2m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{c}_{h}\boldsymbol{r}_{h}^{\top}\right)-\dfrac{\mathbb{E}(U)}{2}\overline{\boldsymbol{c}}_B(\boldsymbol{z})\overline{\boldsymbol{r}}_B(\boldsymbol{z})^\top\\
#'              &\quad+\left(\dfrac{\mathbb{E}(U)}{2m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{r}_{h}\boldsymbol{c}_{h}^{\top}\right)-\dfrac{\mathbb{E}(U)}{2}\overline{\boldsymbol{r}}_B(\boldsymbol{z})\overline{\boldsymbol{c}}_B(\boldsymbol{z})^\top,
#'      \end{aligned}}
#'   where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameters of the latent variables.
#'  \item \code{"General"}: The latent variables do not have any nice properties:
#'      \deqn{\begin{aligned}
#'              \boldsymbol{S}_B(\boldsymbol{z})&=\left(\dfrac{1}{m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{c}_{h}\boldsymbol{c}_{h}^{\top}\right)-\overline{\boldsymbol{c}}_B(\boldsymbol{z})\overline{\boldsymbol{c}}_B(\boldsymbol{z})^\top\\
#'              &\quad+\left(\dfrac{1}{4m}\boldsymbol{\mathfrak{E}}_{UU}\bullet\sum\limits_{h=1}^{n}z_{h}\boldsymbol{r}_{h}\boldsymbol{r}_{h}^{\top}\right)-\dfrac{1}{4}\boldsymbol{\mathfrak{E}}_{UU}\bullet\left[\overline{\boldsymbol{r}}_B(\boldsymbol{z})\overline{\boldsymbol{r}}_B(\boldsymbol{z})^\top\right]\\
#'              &\quad+\left(\dfrac{1}{2m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{c}_{h}\boldsymbol{r}_{h}^{\top}\right)\boldsymbol{\Psi}-\dfrac{1}{2}\overline{\boldsymbol{c}}_B(\boldsymbol{z})\overline{\boldsymbol{r}}_B(\boldsymbol{z})^\top\boldsymbol{\Psi}\\
#'              &\quad+\boldsymbol{\Psi}\left(\dfrac{1}{2m}\sum\limits_{h=1}^{n}z_{h}\boldsymbol{r}_{h}\boldsymbol{c}_{h}^{\top}\right)-\dfrac{1}{2}\boldsymbol{\Psi}\overline{\boldsymbol{r}}_B(\boldsymbol{z})\overline{\boldsymbol{c}}_B(\boldsymbol{z})^\top,
#'      \end{aligned}}
#'    where:
#'      \itemize{
#'          \item \eqn{\boldsymbol{\Psi}=\text{diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt}, 
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'          \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'     }
#' }
#' 
#' @param z A vector of 0 and 1, indicating which observations should be considered for the calculation
#' @param data An \code{\linkS4class{intData}} object containing the macrodata/interval data
#' @return The symbolic covariance matrix
#' @export
#' @examples
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute sample interval-valued covariance matrix using the all the observations
#' z <- rep(1, nrow(credit_card_int))
#' credit_card_cov <- int_cov_z(z, credit_card_int)
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Minimum Covariance Determinant Estimator and Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2604.26769. \url{https://arxiv.org/abs/2604.26769}
int_cov_z <- function(z,data){
    
    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    param <- data@LatentParam
    case <- data@LatentCase
    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    m <- sum(z)

    C_t <- t(C)
    R_t <- t(R)

    zC <- z*C
    zR <- z*R
    sum_c <- C_t%*%zC
    sum_r <- R_t%*%zR
    sum_cr <- C_t%*%zR
    sum_rc <- R_t%*%zC

    zz_t <- z%*%t(z)
    C_z <- C_t%*%zz_t
    R_z <- R_t%*%zz_t

    if (case=="U_id_symmetric"){
        delta <- param[[1]]

        sigma_z <- 1/m*sum_c - (1/m^2)*C_z%*%C +
                + delta/m*sum_r - (delta/m^2)*R_z%*%R
    }else if (case=="U_id") {
        delta <- param[[1]]
        U_mean <- param[[2]]

        sigma_z <- 1/m*sum_c - (1/m^2)*C_z%*%C +
                + delta/m*sum_r - (delta/m^2)*R_z%*%R +
                + U_mean/(2*m)*sum_cr-U_mean/(2*m^2)*C_z%*%R +
                + U_mean/(2*m)*sum_rc-U_mean/(2*m^2)*R_z%*%C
    }else if (case=="General") {
        e_UU <- param[[1]]
        psi <- param[[2]]

        sigma_z <- 1/m*sum_c - (1/m^2)*C_z%*%C +
                + 1/(4*m)*e_UU*sum_r - 1/(4*m^2)*e_UU*R_z%*%R +
                + 1/(2*m)*sum_cr%*%psi-1/(2*m^2)*C_z%*%R%*%psi +
                + 1/(2*m)*psi%*%sum_rc-1/(2*m^2)*psi%*%R_z%*%C
    }
    colnames(sigma_z) <- rownames(sigma_z) <- colnames(data)
    return(sigma_z)
}

#' Safely invert a covariance matrix with Moore-Penrose generalized inverse fallback
#'
#' Computes a numerically stable inverse of a covariance matrix. The function:
#' \enumerate{
#'   \item Attempts standard inversion via \code{solve()}.
#'   \item If the matrix is ill-conditioned, falls back to a Moore-Penrose generalized inverse.
#' }
#'
#' @param cov A numeric covariance matrix.
#' @param verbose Logical; if \code{TRUE} (default), emits warnings when fallback method is used.
#'
#' @return A matrix representing:
#' \itemize{
#'   \item The inverse of \code{cov} if well-conditioned
#'   \item A Moore–Penrose generalized inverse if inversion fails
#' }
#'
#' @details
#' When the covariance matrix is singular or nearly singular, direct inversion
#' may fail or produce unstable results. This function ensures robustness by
#' using Moore-Penrose generalized inverse (via [MASS::ginv()]).
#'
#' The pseudo-inverse effectively ignores directions with negligible variance,
#' which may slightly affect interpretations (e.g., Mahalanobis distances or
#' Shapley values).
#'
#' @examples
#' set.seed(1)
#' 
#' # Example where inversion fails
#' X <- matrix(rnorm(20), ncol = 5)
#' cov_X <- cov(X)
#'
#' #solve(cov_X)  # Standard inversion fails
#' safe_solve_cov(cov_X) # Returns a generalized inverse
#' 
#' # Example where inversion does not fail
#' Y <- cbind(rnorm(20), rnorm(20, mean=1, sd=2))
#' cov_Y <- cov(Y)
#' solve(cov_Y)  # Standard inversion succeeds
#' safe_solve_cov(cov_Y)  # Returns same result
#' @keywords internal
#' @export
safe_solve_cov <- function(cov, verbose = TRUE) {
  
    # Basic checks
    if (!is.matrix(cov)) {
        stop("`cov` must be a matrix.")
    }
    if (nrow(cov) != ncol(cov)) {
        stop("`cov` must be square.")
    }
    
    # Try standard inversion
    inv <- tryCatch(
        solve(cov),
        error = function(e) NULL
    )
    
    # Fallback: generalized inverse (Moore–Penrose)
    if (is.null(inv)) {
        if (verbose) {
            warning("Covariance matrix inversion failed, using Moore-Penrose generalized inverse.")
        }
        
        inv <- MASS::ginv(cov)
    }
    
    return(inv)
}
