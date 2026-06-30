#' Compute Shapley Values for Interval-valued Data
#' 
#' Outlier explanation based on Shapley values for interval-valued data.
#' Decomposes the squared interval-valued Mahalanobis distance into additive outlyingness contributions of
#' the variables.
#'
#' @details
#' The Shapley value decomposes the squared Interval-Mahalanobis distance (see \code{\link{IMah_dist}}) into additive outlyingness contributions of the variables.
#' Let \eqn{\boldsymbol{\mu}_B=(\boldsymbol{\mu}_C^\top,\boldsymbol{\mu}_R^\top)^\top} be the barycenter and \eqn{\boldsymbol{\Sigma}_B} the symbolic covariance matrix (see \code{\link{int_cov}}).
#' The Shapley value of an interval-valued observation \eqn{\boldsymbol{x}=(\boldsymbol{c}^\top,\boldsymbol{r}^\top)^\top}, for the Interval-Mahalanobis distance, is defined according to the \code{LatentCase}:
#' \itemize{
#'  \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{\boldsymbol{\phi}(\boldsymbol{x})=(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right]+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right],}
#'      where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'  \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{\begin{aligned}
#'              \boldsymbol{\phi}(\boldsymbol{x})&=(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right]+\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]\\
#'              &\quad+\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]+\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],
#'     \end{aligned}}
#'      where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#'  \item \code{"General"}: The latent variables do not have any nice properties:
#'      \deqn{\begin{aligned}
#'          \boldsymbol{\phi}(\boldsymbol{x})&=(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right]
#'          +\dfrac{1}{4}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_B^{-1}\right)(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]\\
#'          &\quad+\dfrac{1}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]
#'          +\dfrac{1}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Psi}\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],
#'      \end{aligned}}
#'      where:
#'      \itemize{
#'          \item \eqn{\boldsymbol{\Psi}=\text{Diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt},
#'          \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'          \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'      }
#' }
#' 
#' @param data An \code{\linkS4class{intData}} object containing the interval-valued dataset (macrodata).
#' @param mean_c (Optional) A vector specifying the mean of centers. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function.
#' @param mean_r (Optional) A vector specifying the mean of ranges. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function.
#' @param cov (Optional) A covariance matrix. Defaults to `NULL`, in which case it will be computed 
#' using the \code{\link{IMCD}} function.
#'
#' @return A matrix of Shapley values with row and column names corresponding to the rows and 
#' columns of the input data.
#' 
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Explainable Outlier Detection for Interval-valued Data. 
#' arXiv preprint arXiv:2606.26307. \url{https://arxiv.org/abs/2606.26307}
#' 
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute Shapley values based on IMCD estimates of mean and covariance
#' credit_card_shapley <- int_Shapley(credit_card_int)
#' @export
int_Shapley <- function(data,mean_c=NULL,mean_r=NULL,cov=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase

    if(is.null(mean_c) || is.null(mean_r) || is.null(cov)){
        IMCD_res <- IMCD(data,m=floor(data@NObs*0.75))
        if(is.null(cov)) cov <- IMCD_res$cov_IMCD
        if(is.null(mean_c)) mean_c <- IMCD_res$mean_IMCD_c
        if(is.null(mean_r)) mean_r <- IMCD_res$mean_IMCD_r
    }

    cov_inv <- safe_solve_cov(cov)
    c_0 <- t(C)-mean_c
    r_0 <- t(R)-mean_r

    if (case=="U_id_symmetric"){
        delta <- param[[1]]
        shapley <- c_0*(cov_inv%*%c_0) + delta*r_0*(cov_inv%*%r_0)
    }else if (case=="U_id"){
        delta <- param[[1]]
        U_mean <- param[[2]]
        shapley <- c_0*(cov_inv%*%c_0) + 
                + delta*r_0*(cov_inv%*%r_0) +
                + U_mean/2*c_0*(cov_inv%*%r_0) + 
                + U_mean/2*r_0*(cov_inv%*%c_0)
    }else if (case=="General"){
        e_UU <- param[[1]]
        psi <- param[[2]]
        shapley <- c_0*(cov_inv%*%c_0) + 
                + 1/4*r_0*((e_UU*cov_inv)%*%r_0) +
                + 1/2*c_0*(cov_inv%*%psi%*%r_0) + 
                + 1/2*r_0*(psi%*%cov_inv%*%c_0)
    }
    shapley <- t(shapley)
    rownames(shapley) <- rownames(data)
    colnames(shapley) <- colnames(data)
    return(shapley)
}

#' Compute Shapley Decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) for Interval-valued Data
#' 
#' Decomposes the squared interval-valued Mahalanobis distance of each observation into outlyingness contributions of (Centers, Ranges, and CrossCentersRanges) per variable for interval-valued data.
#' 
#' @details 
#' Let \eqn{\boldsymbol{\mu}_B=(\boldsymbol{\mu}_C^\top,\boldsymbol{\mu}_R^\top)^\top} be the barycenter and \eqn{\boldsymbol{\Sigma}_B} the symbolic covariance matrix (see \code{\link{int_cov}}).
#' Based on the Shapley value (see \code{\link{int_Shapley}}), we can further decompose the Interval-Mahalanobis distance of an interval-valued observation \eqn{\boldsymbol{x}=(\boldsymbol{c}^\top,\boldsymbol{r}^\top)^\top} into contributions of the centers, ranges and cross-centers-ranges of the variables. The decomposition is defined according to the \code{LatentCase}:
#' \itemize{
#'      \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'            \itemize{
#'                 \item Centers contribution: \deqn{(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],}
#'                 \item Ranges contribution: \deqn{\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right],}
#'             }
#'         where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'       \item \code{"U_id"}: The latent variables are identically distributed:
#'             \itemize{
#'                 \item Centers contribution: \deqn{(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],}
#'                 \item Ranges contribution: \deqn{\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right],}
#'                 \item CrossCentersRanges contribution: 
#'                      \deqn{\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]+\dfrac{\mathbb{E}(U)}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],}
#'             }
#'         where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#'       \item \code{"General"}: The latent variables do not have any nice properties:
#'             \itemize{
#'                 \item Centers contribution: \deqn{(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],}
#'                 \item Ranges contribution: \deqn{\dfrac{1}{4}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\left(\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_B^{-1}\right)(\boldsymbol{r}-\boldsymbol{\mu}_R)\right],}
#'                 \item CrossCentersRanges contribution: 
#'                      \deqn{\dfrac{1}{2}(\boldsymbol{c}-\boldsymbol{\mu}_C)\bullet\left[\boldsymbol{\Sigma}_B^{-1}\boldsymbol{\Psi}(\boldsymbol{r}-\boldsymbol{\mu}_R)\right]+\dfrac{1}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)\bullet\left[\boldsymbol{\Psi}\boldsymbol{\Sigma}_B^{-1}(\boldsymbol{c}-\boldsymbol{\mu}_C)\right],}
#'             }
#'          where:
#'          \itemize{
#'              \item \eqn{\boldsymbol{\Psi}=\text{Diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt},
#'              \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'              \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#'          } 
#' }
#'    
#' @param data An \code{\linkS4class{intData}} object containing the interval-valued dataset (macrodata).
#' @param mean_c (Optional) A vector specifying the mean of centers. Defaults to `NULL`, in which case it will be computed using the \code{\link{IMCD}} function.
#' @param mean_r (Optional) A vector specifying the mean of ranges. Defaults to `NULL`, in which case it will be computed using the \code{\link{IMCD}} function.
#' @param cov (Optional) A covariance matrix. Defaults to `NULL`, in which case it will be computed using the \code{\link{IMCD}} function.
#' 
#' @return A list containing the matrix of Shapley value decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) per variable for each observation.
#' 
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Explainable Outlier Detection for Interval-valued Data.
#' arXiv preprint arXiv:2606.26307. \url{https://arxiv.org/abs/2606.26307}
#' 
#' @export
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute Shapley decomposition into contributions of (Centers, Ranges, and CrossCentersRanges) 
#' # based on IMCD estimates of mean and covariance
#' credit_card_shap_decomp <- int_Shapley_decomp(credit_card_int)
int_Shapley_decomp <- function(data,mean_c=NULL,mean_r=NULL,cov=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase

    if(is.null(mean_c) || is.null(mean_r) || is.null(cov)){
        IMCD_res <- IMCD(data,m=floor(data@NObs*0.75))
        if(is.null(cov)) cov <- IMCD_res$cov_IMCD
        if(is.null(mean_c)) mean_c <- IMCD_res$mean_IMCD_c
        if(is.null(mean_r)) mean_r <- IMCD_res$mean_IMCD_r
    }

    cov_inv <- safe_solve_cov(cov)
    shapley_decomp <- list()

    for(i in 1:data@NObs){
        c_0 <- C[i,]-mean_c
        r_0 <- R[i,]-mean_r

        if (case=="U_id_symmetric"){
            delta <- param[[1]]
            shapley <- cbind(c_0*(cov_inv%*%c_0), delta*r_0*(cov_inv%*%r_0))
            colnames(shapley) <- c("Centers","Ranges")
        }else if (case=="U_id"){
            delta <- param[[1]]
            U_mean <- param[[2]]
            shapley <- cbind(c_0*(cov_inv%*%c_0), 
                    delta*r_0*(cov_inv%*%r_0),
                    U_mean/2*c_0*(cov_inv%*%r_0)+U_mean/2*r_0*(cov_inv%*%c_0))
            colnames(shapley) <- c("Centers","Ranges","CentersRanges")
        }else if (case=="General"){
            e_UU <- param[[1]]
            psi <- param[[2]]
            shapley <- cbind(c_0*(cov_inv%*%c_0),
                    1/4*r_0*((e_UU*cov_inv)%*%r_0),
                    1/2*c_0*(cov_inv%*%psi%*%r_0)+1/2*r_0*(psi%*%cov_inv%*%c_0))
            colnames(shapley) <- c("Centers","Ranges","CentersRanges")
        }
        rownames(shapley) <- colnames(data)
        shapley_decomp[[rownames(data)[i]]] <- t(shapley)
    }
    return(shapley_decomp)
}

#' Compute Shapley interaction indices for Interval-valued Data
#' 
#' Obtains a \eqn{p \times p} matrix containing pairwise outlyingness scores based on Shapley interaction indices for each observation.
#' Decomposes the squared interval-valued Mahalanobis distance of each observation into outlyingness contributions of pairs of variables.
#'
#' @details
#' Let \eqn{\boldsymbol{\mu}_B=(\boldsymbol{\mu}_C^\top,\boldsymbol{\mu}_R^\top)^\top} be the barycenter and \eqn{\boldsymbol{\Sigma}_B} the symbolic covariance matrix (see \code{\link{int_cov}}).
#' Let also \eqn{\boldsymbol{\phi}(\boldsymbol{x})} be the Shapley value of \eqn{\boldsymbol{x}} (see \code{\link{int_Shapley}}) and \eqn{\mathrm{diag}(\boldsymbol{v})} be the diagonal matrix whose main diagonal is the vector \eqn{\boldsymbol{v}}.
#' The Shapley interaction index of an interval-valued observation \eqn{\boldsymbol{x}=(\boldsymbol{c}^\top,\boldsymbol{r}^\top)^\top}, for the Interval-Mahalanobis distance, is defined according to the \code{LatentCase}:
#' \itemize{
#'  \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric:
#'      \deqn{\boldsymbol{\Phi}(\boldsymbol{x})=2(\boldsymbol{c}-\boldsymbol{\mu}_C)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\bullet\boldsymbol{\Sigma}_B^{-1} + 2\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)(\boldsymbol{r}-\boldsymbol{\mu}_R)^\top\bullet\boldsymbol{\Sigma}_B^{-1}-\mathrm{diag}\left(\boldsymbol{\phi}(\boldsymbol{x})\right),}
#'      where \eqn{\delta=\mathbb{E}(U^2)/4} is the parameter of the latent variables.
#'  \item \code{"U_id"}: The latent variables are identically distributed:
#'     \deqn{\begin{aligned}
#'             \boldsymbol{\Phi}(\boldsymbol{x})&=2(\boldsymbol{c}-\boldsymbol{\mu}_C)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\bullet\boldsymbol{\Sigma}_B^{-1} + 2\delta(\boldsymbol{r}-\boldsymbol{\mu}_R)(\boldsymbol{r}-\boldsymbol{\mu}_R)^\top\bullet\boldsymbol{\Sigma}_B^{-1}\\
#'            &\quad+\mathbb{E}(U)(\boldsymbol{c}-\boldsymbol{\mu}_C)(\boldsymbol{r}-\boldsymbol{\mu}_R)^\top\bullet\boldsymbol{\Psi} + \mathbb{E}(U)(\boldsymbol{r}-\boldsymbol{\mu}_R)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\bullet\boldsymbol{\Sigma}_B^{-1}-\mathrm{diag}\left(\boldsymbol{\phi}(\boldsymbol{x})\right),
#'    \end{aligned}}
#'    where \eqn{\delta=\mathbb{E}(U^2)/4} and \eqn{\mathbb{E}(U)} are the parameter of the latent variables.
#' \item \code{"General"}: The latent variables do not have any nice properties:
#'    \deqn{\begin{aligned}
#'       \boldsymbol{\Phi}(\boldsymbol{x})&=2(\boldsymbol{c}-\boldsymbol{\mu}_C)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\bullet\boldsymbol{\Sigma}_B^{-1} + \dfrac{1}{2}(\boldsymbol{r}-\boldsymbol{\mu}_R)(\boldsymbol{r}-\boldsymbol{\mu}_R)^\top\bullet\boldsymbol{\mathfrak{E}}_{UU}\bullet\boldsymbol{\Sigma}_B^{-1}\\
#'      &\quad+(\boldsymbol{c}-\boldsymbol{\mu}_C)(\boldsymbol{r}-\boldsymbol{\mu}_R)^\top\bullet\boldsymbol{\Sigma}_B^{-1}\boldsymbol{\Psi} + (\boldsymbol{r}-\boldsymbol{\mu}_R)(\boldsymbol{c}-\boldsymbol{\mu}_C)^\top\bullet\boldsymbol{\Psi}\boldsymbol{\Sigma}_B^{-1}-\mathrm{diag}\left(\boldsymbol{\phi}(\boldsymbol{x})\right),
#'   \end{aligned}}
#'   where:
#'    \item \eqn{\boldsymbol{\Psi}=\text{Diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))},
#'    \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt},
#'    \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p},
#'    \item \eqn{\bullet} denotes the Schur (or entrywise) product of matrices.
#' }
#' @param data An \code{\linkS4class{intData}} object containing the interval-valued dataset (macrodata).
#' @param mean_c (Optional) A vector specifying the mean of centers. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function.
#' @param mean_r (Optional) A vector specifying the mean of ranges. Defaults to `NULL`, in which case 
#' it will be computed using the \code{\link{IMCD}} function.
#' @param cov (Optional) A covariance matrix. Defaults to `NULL`, in which case it will be computed 
#' using the \code{\link{IMCD}} function.
#'
#' @return A list containing the matrix of Shapley interaction indices for each observation.
#' 
#' @references Loureiro, C. P., Oliveira, M. R., Brito, P., & Oliveira, L. (2026). 
#' Explainable Outlier Detection for Interval-valued Data.
#' arXiv preprint arXiv:2606.26307. \url{https://arxiv.org/abs/2606.26307}
#' 
#' @examples
#' # Create intData object
#' data(creditcard)
#' credit_card_int <- creditcard$intData
#' 
#' # Compute Shapley interaction indices based on the mean and covariance matrix estimated by IMCD
#' credit_card_shap_inter <- int_Shapley_interaction(credit_card_int)
#' @export
int_Shapley_interaction <- function(data,mean_c=NULL,mean_r=NULL,cov=NULL){

    if(!inherits(data,"intData")) stop("Argument data is not an object of class intData\n")

    C <- as.matrix(data@Centers)
    R <- as.matrix(data@Ranges)
    param <- data@LatentParam
    case <- data@LatentCase

    if(is.null(mean_c) || is.null(mean_r) || is.null(cov)){
        IMCD_res <- IMCD(data,m=floor(data@NObs*0.75))
        if(is.null(cov)) cov <- IMCD_res$cov_IMCD
        if(is.null(mean_c)) mean_c <- IMCD_res$mean_IMCD_c
        if(is.null(mean_r)) mean_r <- IMCD_res$mean_IMCD_r
    }

    cov_inv <- safe_solve_cov(cov)
    shapley_interaction <- list()

    for(i in 1:data@NObs){
        c_0 <- C[i,]-mean_c
        r_0 <- R[i,]-mean_r
        if (case=="U_id_symmetric"){
            delta <- param[[1]]
            inter_shapley <- 2*tcrossprod(c_0)*cov_inv + 2*delta*tcrossprod(r_0)*cov_inv
        }else if (case=="U_id"){
            delta <- param[[1]]
            U_mean <- param[[2]]
            inter_shapley <- 2*tcrossprod(c_0)*cov_inv + 2*delta*tcrossprod(r_0)*cov_inv + 
                    + U_mean*tcrossprod(c_0,r_0)*cov_inv + U_mean*tcrossprod(r_0,c_0)*cov_inv
        }else if (case=="General"){
            e_UU <- param[[1]]
            psi <- param[[2]]
            inter_shapley <- 2*tcrossprod(c_0)*cov_inv + 1/2*tcrossprod(r_0)*e_UU*cov_inv +
                    + tcrossprod(c_0,r_0)*(cov_inv%*%psi) + tcrossprod(r_0,c_0)*(psi%*%cov_inv)
        }
        inter_shapley <- inter_shapley - diag(rowSums(inter_shapley/2))
        rownames(inter_shapley) <- colnames(inter_shapley) <- colnames(data)
        shapley_interaction[[rownames(data)[i]]] <- inter_shapley
    }
    return(shapley_interaction)
}
