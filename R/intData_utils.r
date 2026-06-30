#' Compute Latent Variables
#' 
#' Obtain the latent variables inherent to the macrodata.
#' 
#' @details
#' The latent variables, \eqn{U_{j}}, are defined according to the following model:
#' 
#' Let \eqn{X_j=(C_j,R_j)^\top=\left[C_j-\dfrac{R_j}{2}, C_j+\dfrac{R_j}{2}\right]} represent the \strong{macrodata} and
#' \deqn{V_{j}=C_j+U_{j}\dfrac{R_j}{2},\quad j=1,\dots,p,}
#' the \strong{microdata} with \eqn{U_{j}} being random variables with support on \eqn{[-1,1]}, uncorrelated with \eqn{(C_j,R_j)}.
#' 
#' @param microdata A matrix containing the microdata.
#' @param macrodata A data frame, matrix or \code{\linkS4class{intData}} object containing the macrodata/interval data.
#' @param agrby A factor used to specify the grouping of the microdata.
#' @param agrlevels The categories/levels on which the microdata was aggregated.
#' @param Seq Format of macrodata if it is a data frame or matrix. Available options are:
#' \itemize{
#'   \item \code{"AllLb_AllUb"}: All lower bounds followed by all upper bounds, in the same variable order.
#'   \item \code{"AllCen_AllRng"}: All Centers followed by all Ranges, in the same variable order.
#'   \item \code{"LbUb_VarbyVar"}: Lower bounds followed by upper bounds, variable by variable.
#'   \item \code{"CenRng_VarbyVar"}: Centers followed by Ranges, variable by variable.
#' }
#' @return A matrix with the same size as the microdata.
#' @references Oliveira, M.R., Azeitona, M., Pacheco, A., Valadas, R.. Association measures for interval variables. Advances in Data Analysis and Classification 16, 491–520 (2022). \doi{10.1007/s11634-021-00445-8}
#' @export
#' @examples
#' data(creditcard)
#' CreditCard_min_max <- creditcard$min_max
#' CreditCard_microdata <- creditcard$microdata
#' 
#' # Define grouping variable for microdata aggregation
#' credit_agrby <- paste(CreditCard_microdata$Name, CreditCard_microdata$Month, sep = "_")
#' 
#' # Obtain latent variables inherent to the macrodata (standardized to [-1,1])
#' credit_card_U <- get_latent_var(microdata = CreditCard_microdata[,3:7], 
#'                                 macrodata = CreditCard_min_max, 
#'                                 agrby = credit_agrby, 
#'                                 agrlevels = row.names(CreditCard_min_max), 
#'                                 Seq = "LbUb_VarbyVar")
get_latent_var <- function(microdata,
                            macrodata,
                            agrby,
                            agrlevels,
                            Seq=c("AllLb_AllUb","AllCen_AllRng","LbUb_VarbyVar","CenRng_VarbyVar")){
    Seq <- match.arg(Seq)

    if(inherits(macrodata, "intData")){
        C <- macrodata@Centers
        R <- macrodata@Ranges
    }else if(is.data.frame(macrodata) || is.matrix(macrodata)){
        p <- ncol(macrodata)
        q <- p/2
        if (Seq == "LbUb_VarbyVar") {Lbnd <- macrodata[,2*(0:(q-1))+1]; Ubnd <- macrodata[,2*(1:q)]; C <- (Lbnd+Ubnd)/2; R <- (Ubnd-Lbnd)}
        if (Seq == "AllLb_AllUb")   {Lbnd <- macrodata[,1:q] ; Ubnd <- macrodata[,(q+1):p]; C <- (Lbnd+Ubnd)/2; R <- (Ubnd-Lbnd)}
        if (Seq == "CenRng_VarbyVar") {C <- macrodata[,2*(0:(q-1))+1]; R <- macrodata[,2*(1:q)]}
        if (Seq == "AllCen_AllRng")   {C <- macrodata[,1:q]; R <- macrodata[,(q+1):p]}
    }else stop("macrodata must be a data.frame, matrix, or intData object.")

    if (!is.data.frame(C)) C <- as.data.frame(C)
    if (!is.data.frame(R)) R <- as.data.frame(R)

    U <- microdata
    for (i in seq_len(nrow(U))){
        for (j in seq_len(ncol(U))){
            c_ij <- C[agrlevels==agrby[i],j]
            r_ij <- R[agrlevels==agrby[i],j]
            U[i,j] <- 2*(microdata[i,j]-c_ij)/r_ij
        }
    }

    U[U==1|U==-1] <- NA
    return(U)
}

#' Compute Latent Variables Parameters
#' 
#' Obtain the parameters of the latent variables inherent to the macrodata.
#' 
#' @details 
#' The parameters of the latent variables inherent to the macrodata are defined according to the \code{LatentCase}:
#' \itemize{
#'  \item \code{"U_id_symmetric"}: The latent variables are identically distributed and symmetric, so its parameters are:
#'      \itemize{
#'          \item  \eqn{\delta=\mathbb{E}(U^2)/4}
#'      }
#'  \item \code{"U_id"}: The latent variables are identically distributed, so its parameters are:
#'      \itemize{
#'         \item  \eqn{\delta=\mathbb{E}(U^2)/4}
#'         \item  \eqn{\mathbb{E}(U)}
#'      }
#'  \item \code{"General"}: The latent variables do not have any nice properties, so its parameters are:
#'      \itemize{
#'         \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt}, and \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p}
#'         \item \eqn{\boldsymbol{\Psi}=\text{Diag}(\mathbb{E}(U_1),\dots,\mathbb{E}(U_p))}
#'      }
#' }
#' @param LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"}, \code{"Triang"}, \code{"TNorm"}, \code{"InvTri"}, \code{"Beta"}, \code{"KDE"}, \code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#'  The default is \code{"Unif"} if \code{LatentCase="U_id_symmetric"} or if \code{Umicro} is not provided, and \code{"KDE"} if \code{LatentCase="General"}.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param Umicro Latent microdata observations.  
#'   Needed if \code{estimate.DistParam} is \code{TRUE} or \code{LatentDist} is \code{"KDE"}.
#' @param p Number of variables.
#' @param estimate.DistParam Logical parameter indicating if estimation of the parameters of the latent distributions should be performed. Can only be set to TRUE if \code{LatentCase="General"}.
#' The default is \code{FALSE}.
#' @return A list with the parameters of the latent variables.
#' @export
#' @examples
#' data(creditcard)
#' CreditCard_min_max <- creditcard$min_max
#' CreditCard_microdata <- creditcard$microdata
#' 
#' # Define grouping variable for microdata aggregation
#' credit_agrby <- paste(CreditCard_microdata$Name, CreditCard_microdata$Month, sep = "_")
#' 
#' # Obtain latent variables inherent to the macrodata (standardized to [-1,1])
#' credit_card_U <- get_latent_var(microdata = CreditCard_microdata[,3:7], 
#'                                 macrodata = CreditCard_min_max, 
#'                                 agrby = credit_agrby, 
#'                                 agrlevels = row.names(CreditCard_min_max), 
#'                                 Seq = "LbUb_VarbyVar")
#' 
#' # Obtain parameters of the latent variables
#' credit_card_param <- get_latent_param(LatentCase = "General",
#'                                       LatentDist = "KDE",
#'                                       Umicro = credit_card_U)
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
get_latent_param <- function(LatentCase=c("U_id_symmetric","U_id","General"),
                             LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                             TriangParam=0,
                             BetaParam.a=1,
                             BetaParam.b=1,
                             Umicro=NULL,
                             p=NULL,
                             estimate.DistParam=FALSE){

    case <- match.arg(LatentCase)
    
    if(case!="General"&&length(unique(TriangParam)) > 1) stop("Error: For different TriangParam for each variable, LatentCase must be 'General'.")
    if(case!="General"&&length(unique(BetaParam.a)) > 1) stop("Error: For different BetaParam.a for each variable, LatentCase must be 'General'.")

    if (!identical(LatentDist, c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"))){dist <- match.arg(LatentDist,several.ok = TRUE)}
    else {
        if(case=="U_id") stop("Error: a distribution type must be provided for `LatentCase='U_id'`.")
        else if(case=="U_id_symmetric"||is.null(Umicro)) dist <- "Unif"
        else dist <- "KDE"
    }

    if(!is.null(Umicro)) p <- ncol(Umicro)
    if(case!="U_id_symmetric" && is.null(p) && is.null(Umicro)) stop("The number of variables or the microdata must be provided.")
    if(any(dist == "KDE") && is.null(Umicro)) stop("For distribution type KDE, the microdata must be provided.")
    if(length(dist) > 1 && length(dist) != p) stop("Error: a distribution must be provided for each variable or a single distribution must be provided to use for all variables.")
    if (length(dist)==1 && (dist=="KDE" || length(TriangParam)==p || length(BetaParam.a)==p)) dist <- rep(dist,p)

    if(case=="General"&&estimate.DistParam&&is.null(Umicro)) stop("To estimate the distribution parameters, the microdata must be provided.")
    if(case=="General"&&estimate.DistParam){
        TriangParam <- BetaParam.a <- BetaParam.b <- rep(NA,p)
        for(i in 1:p){
            if(dist[i]=="Triang") TriangParam[i] <- 3*mean(Umicro[,i],na.rm=TRUE)
            else if(dist[i]=="Beta") {
                fit <- MASS::fitdistr((Umicro[,i]+1)/2, "beta", start = list(shape1 = 1, shape2 = 1))
                BetaParam.a[i] <- fit$estimate["shape1"]
                BetaParam.b[i] <- fit$estimate["shape2"]
            }
        }
    }else if(case!="General"&&estimate.DistParam) stop("To estimate the distribution parameters, set the LatentCase to 'General'.")

    if(case=="U_id_symmetric"){
        if(dist=="InvTri")      {delta <- 1/8}
        else if(dist=="Unif")   {delta <- 1/12}
        else if(dist=="Triang") {delta <- 1/24}
        else if(dist=="TNorm")  {delta <- 1/36-dnorm(3)/(12*pnorm(3)-6)}
        else if(dist=="Degenerated"){
            message("The latent variable, U, is a degenerated random variable, i.e., P(U=0)=1 and the data is considered as non-symbolic")
            delta <- 0
        }else{
            stop("Not an admissable distribution type for this latent case.")
        }
        param <- list(delta)
    }else if(case=="U_id"){
        delta <- meanU2(dist,TriangParam,BetaParam.a,BetaParam.b,Umicro,p)/4
        U_mean <- meanU(dist,TriangParam,BetaParam.a,BetaParam.b,Umicro,p)
        param <- list(delta,U_mean)
    }else{
        e_UU <- cal.E.UU(dist,TriangParam,BetaParam.a,BetaParam.b,Umicro,p)
        psi <- meanU(dist,TriangParam,BetaParam.a,BetaParam.b,Umicro,p)
        param <- list(e_UU,psi)
    }
    if (length(unique(dist)) == 1) dist <- dist[1] 
    return(list(LatentParam=param,TriangParam=TriangParam,BetaParam.a=BetaParam.a,BetaParam.b=BetaParam.b,LatentCase=case,LatentDist=dist))
}

#' Compute Mean Latent Variables
#' 
#' Obtain the mean of the latent variables inherent to the macrodata.
#' 
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"}, \code{"Triang"}, \code{"TNorm"}, \code{"InvTri"}, \code{"Beta"}, \code{"KDE"}, \code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param Umicro Latent microdata observations. Needed if \code{LatentDist="KDE"}.
#' @param p Number of variables.
#' @return Either a diagonal matrix with the mean of each variable or a value if the variables are identically distributed.
#' @keywords internal
meanU <- function(LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                    TriangParam=0,
                    BetaParam.a=1,
                    BetaParam.b=1,
                    Umicro=NULL,
                    p=NULL){
    dist <- match.arg(LatentDist,several.ok = TRUE)

    if(!is.null(Umicro)) p <- ncol(Umicro)
    if (length(dist)==1 && (dist=="KDE" || length(TriangParam)==p || length(BetaParam.a)==p)) dist <- rep(dist,p)

    if(length(dist)==1){
        if(dist=="Unif") mU <- 0
        else if(dist=="Beta") mU <- 2*BetaParam.a/(BetaParam.a+BetaParam.b)-1
        else if(dist=="Triang") mU <- TriangParam/3
        else stop("Not an admissable distribution type")
    }else if(length(dist)==p){
        mU <- rep(NA,p)
        for(i in 1:p){
            if(dist[i] == "KDE") {mU[i] <- mean(Umicro[,i],na.rm=TRUE)}
            else if(dist[i] == "Unif") {mU[i] <- 0}
            else if(dist[i] == "Beta") {mU[i] <- 2*BetaParam.a[i]/(BetaParam.a[i]+BetaParam.b[i])-1}
            else if(dist[i] == "Triang") {mU[i] <- TriangParam[i]/3}
            else stop("Not an admissable distribution type")
        }
        mU <- diag(mU)
    }else{
       stop("Error: a distribution must be provided for each variable.")
    }
    return(mU)
}

#' Compute Mean Square Latent Variables
#' 
#' Obtain the mean of the square of the latent variables inherent to the macrodata.
#' 
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"}, \code{"Triang"}, \code{"TNorm"}, \code{"InvTri"}, \code{"Beta"}, \code{"KDE"}, \code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param Umicro Latent microdata observations. Needed if \code{LatentDist="KDE"}.
#' @param p Number of variables.
#' @return Either a diagonal matrix with the mean of the square of each variable or a value if the variables are identically distributed.
#' @keywords internal
meanU2 <- function(LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                    TriangParam=0,
                    BetaParam.a=1,
                    BetaParam.b=1,
                    Umicro=NULL,
                    p=NULL){
    dist <- match.arg(LatentDist,several.ok = TRUE)

    if(!is.null(Umicro)) p <- ncol(Umicro)
    if (length(dist)==1 && (dist=="KDE" || length(TriangParam)==p || length(BetaParam.a)==p)) dist <- rep(dist,p)

    if(length(dist)==1){
        if(dist=="Unif") mU2 <- 1/3
        else if(dist=="Beta") mU2 <- 4*BetaParam.a*BetaParam.b/( (BetaParam.a+BetaParam.b)^2*(BetaParam.a+BetaParam.b-1) ) + abs(BetaParam.a-BetaParam.b)/(BetaParam.a+BetaParam.b)
        else if(dist=="Triang") mU2 <- (1+TriangParam^2)/6
        else stop("Not an admissable distribution type")
    }else if(length(dist)==p){
        mU2 <- rep(NA,p)
        for(i in 1:p){
            if(dist[i] == "KDE") mU2[i] <- mean(Umicro[,i]^2,na.rm=TRUE)
            else if(dist[i] == "Unif") mU2[i] <- 1/3
            else if(dist[i] == "Beta") mU2[i] <- 4*BetaParam.a[i]*BetaParam.b[i]/( (BetaParam.a[i]+BetaParam.b[i])^2*(BetaParam.a[i]+BetaParam.b[i]-1) ) + abs(BetaParam.a[i]-BetaParam.b[i])/(BetaParam.a[i]+BetaParam.b[i])
            else if(dist[i] == "Triang") mU2[i] <- (1+TriangParam[i]^2)/6
            else stop("Not an admissable distribution type")
        }
        mU2 <- diag(mU2)
    }else{
       stop("Error: a distribution must be provided for each variable.")
    }
    return(mU2)
}

#' Compute Cal.E Latent Variables
#' 
#' Computes \eqn{\boldsymbol{\mathfrak{E}}_{UU}} for the latent variables inherent to the macrodata.
#' 
#' @details
#' The matrix \eqn{\boldsymbol{\mathfrak{E}}_{UU}} is defined as follows:
#' \itemize{
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{j\ell}=\mathcal{E}(U_j,U_\ell)}, \eqn{j\neq \ell}, with \eqn{\mathcal{E}(U_j,U_\ell)=\int_0^1 F_{U_j}^{-1}(t) F_{U_\ell}^{-1}(t) \, dt}
#'      \item \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{jj}=\mathbb{E}(U_j^2)}, \eqn{j,\ell=1,\dots,p}.
#' }
#' 
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"}, \code{"Triang"}, \code{"TNorm"}, \code{"InvTri"}, \code{"Beta"}, \code{"KDE"}, \code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param Umicro Latent microdata observations. Needed if \code{LatentDist="KDE"}.
#' @param p Number of variables.
#' @return A \eqn{p\times p} matrix.
#' @keywords internal
cal.E.UU <- function(LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                        TriangParam=0,
                        BetaParam.a=1,
                        BetaParam.b=1,
                        Umicro=NULL,
                        p=NULL){
    dist <- match.arg(LatentDist,several.ok = TRUE)
    if(!is.null(Umicro)) p <- ncol(Umicro)

    if(length(dist)==1) dist <- rep(dist,p)
    if(length(TriangParam)==1)  TriangParam <- rep(TriangParam,p)
    if(length(BetaParam.a)==1)  BetaParam.a <- rep(BetaParam.a,p)
    if(length(BetaParam.b)==1)  BetaParam.b <- rep(BetaParam.b,p)

    cal.E <- meanU2(dist,TriangParam,BetaParam.a,BetaParam.b,Umicro,p)
    
    for(i in 1:(p-1)){
        x <- Umicro[!is.na(Umicro[,i]),i]
        for(j in (i+1):p){
            y <- Umicro[!is.na(Umicro[,j]),j]
            if((dist[i] == "KDE") && (dist[j] == "KDE")){
                cal.E[i,j] <- cal.E[j,i] <- CalE.kde.kde(x,y)
            }
            else if( (dist[i] == "Beta") && (dist[j] == "KDE")){
                cal.E[i,j] <- cal.E[j,i] <- CalE.beta.kde(y,a1=BetaParam.a[i],b1=BetaParam.b[i])
            }
            else if( (dist[i] == "KDE") && (dist[j] == "Beta")){
                cal.E[i,j] <- cal.E[j,i] <- CalE.beta.kde(x,a1=BetaParam.a[j],b1=BetaParam.b[j])
            }
            else if( (dist[i] == "Beta") && (dist[j] == "Beta")){
                cal.E[i,j] <- cal.E[j,i] <- CalE.beta.beta(a1=BetaParam.a[i],b1=BetaParam.b[i],a2=BetaParam.a[j],b2=BetaParam.b[j])
            }
            else if( (dist[i] == "Triang") && (dist[j] == "Triang")){
                cal.E[i,j] <- cal.E[j,i] <- CalE.triang.triang(mo1=TriangParam[i],mo2=TriangParam[j])
            }
            else{
                stop("Not an admissable combination of two microdata distribution types")
            }
        }
    }
    return(cal.E)
}

#' Computes \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)} for the latent variables inherent to the macrodata, where they follow a Beta distribution.
#' @param a1 Parameter alpha of the first Beta distribution.
#' @param a2 Parameter alpha of the second Beta distribution.
#' @param b1 Parameter beta of the first Beta distribution.
#' @param b2 Parameter beta of the second Beta distribution.
#' @return  Value
#' @keywords internal
CalE.beta.beta <- function(a1,b1,a2,b2){
    integrandBeta <- function(x,a1,b1,a2,b2) {qbeta(x,a1,b1)*qbeta(x,a2,b2)}

    calE.aux <- integrate(integrandBeta, lower=0, upper=1, a1=a1, b1=b1, a2=a2, b2=b2)$value
    calE12 <- 1 - 2*a1/(a1+b1) - 2*a2/(a2+b2) + 4*calE.aux

    return(calE12)
}

#' Computes \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)} for the latent variables inherent to the macrodata, where U_1 follows a Beta(a_1,b_1) and the PDF of U_2 is estimated by a KDE.
#' @param micro Latent microdata observations.
#' @param a1 Parameter alpha of the Beta distribution.
#' @param b1 Parameter beta of the Beta distribution.
#' @return  Value
#' @keywords internal
CalE.beta.kde <- function(micro,a1,b1){
    fit3 <- kde1d::kde1d(micro) # estimate density
    integrandBetaBeta.kde <- function(x,fit3,a=a1,b=b1) {qbeta(x,a,b)*kde1d::qkde1d(x, fit3)}

    calE.aux <- integrate(integrandBetaBeta.kde, lower = 0, upper = 1,fit3)$value
    calE13 <- 2*calE.aux - mean(micro)  

    return(calE13)
}

#' Computes \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)} for the latent variables inherent to the macrodata, where the PDF is estimated by a KDE.
#' @param micro1 Latent microdata observations of the first latent variable.
#' @param micro2 Latent microdata observations of the second latent variable.
#' @return  Value
#' @keywords internal
CalE.kde.kde <- function(micro1,micro2){
    fit3 <- kde1d::kde1d(micro1) # estimate density
    fit4 <- kde1d::kde1d(micro2) # estimate density
    integrand.kde.kde <- function(x,fit3=fit3,fit4=fit4) {kde1d::qkde1d(x, fit3)*kde1d::qkde1d(x, fit4)}

    calE.aux <- integrate(integrand.kde.kde, lower=0, upper=1, fit3, fit4)$value

    return(calE.aux)
}

#' Computes \eqn{[\boldsymbol{\mathfrak{E}}_{UU}]_{ij}=\mathcal{E}(U_i,U_j)} for the latent variables inherent to the macrodata, where they follow a Triangular distribution.
#' @param mo1 Mode of the triangular distribution of the first latent variable.
#' @param mo2 Mode of the triangular distribution of the second latent variable.
#' @return  Value
#' @keywords internal
CalE.triang.triang <- function(mo1=0,mo2=0){
    md1 <- min(mo1,mo2)
    md2 <- max(mo1,mo2)

    if(md1==0 && md2==0)  {calE.aux <- 1/6}
    else{
        integrand1 <- function(x,m1=md1,m2=md2) {(-1+sqrt(2*x*(m1+1))) * (-1+sqrt(2*(1+m2)*x))}
        calE.aux <-  integrate(integrand1, lower=0, upper=(md1+1)/2, m1=md1, m2=md2)$value
    }
    return(calE.aux)
}