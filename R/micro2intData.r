#' Aggregate Microdata into Interval-Valued Data
#'
#' Aggregates microdata from a data frame into interval-valued data using various criteria and latent distribution settings.
#'
#' @param MicDtDF A data frame containing the microdata. All columns should be numeric.
#' @param agrby A factor used to specify the grouping of the microdata for aggregation.
#' @param agrcrt A string or numeric vector of length 2 specifying the aggregation criterion. The default is \code{"minmax"}, which takes the minimum and maximum values for each variable. If a numeric vector is provided, it should specify the lower and upper percentiles for aggregation (e.g., \code{c(0.05, 0.95)}).
#' @param LatentParam Optional latent parameter used for certain types of latent distributions.
#' @param LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"},\code{"Triang"},\code{"TNorm"},\code{"InvTri"},\code{"Beta"},\code{"KDE"},\code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' The default is \code{"KDE"} if \code{LatentCase="General"}.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param estimate.DistParam Logical parameter indicating if estimation of the parameters of the latent distributions should be performed. Can only be set to TRUE if \code{LatentCase="General"}.
#' The default is \code{FALSE}.
#'
#' @return An \linkS4class{intData} object containing the aggregated interval-valued data, or \code{NULL} if all units lead to degenerate intervals.
#'
#' @importFrom assertthat is.number
#' @importFrom stats quantile
#' 
#' @details
#' This function processes a data frame of microdata and aggregates it into interval-valued data according to the specified grouping factor and aggregation criteria. 
#' It can handle different latent distribution cases and parameter settings.
#' 
#' If some rows contain invalid (non-finite or missing) values, those rows are removed before aggregation. If all rows in the resulting interval-valued data are degenerate (i.e., the lower bound equals the upper bound), the function will return \code{NULL}.
#' 
#' @references Adapted from package \code{MAINT.Data} (\url{https://cran.r-project.org/package=MAINT.Data}).
#' 
#' @examples
#' data(creditcard)
#' CreditCard_microdata <- creditcard$microdata
#' credit_agrby<-factor(paste(CreditCard_microdata$Name,CreditCard_microdata$Month,sep = "_"))
#' credit_agr<-micro2intData(CreditCard_microdata[,3:7],credit_agrby,LatentCase = "General")
#' 
#' @export
micro2intData <- function(MicDtDF,
                          agrby,
                          agrcrt="minmax",
                          LatentParam=NULL,
                          LatentCase=c("U_id_symmetric","U_id","General"),
                          LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                          TriangParam=0,
                          BetaParam.a=1,
                          BetaParam.b=1,
                          estimate.DistParam=FALSE){
  mcall <- match.call()$MicDtDF
  if (length(mcall) > 1) mcall <- "microdata"
  if (!(is.data.frame(MicDtDF))) stop("First argument of AgMicroData must be a data frame.\n")
  if (!is.data.frame(MicDtDF)) MicDtDF <- as.data.frame(MicDtDF)
  if (any(!sapply(seq_len(ncol(MicDtDF)),function(ind) is.numeric(MicDtDF[,ind])))){  
    stop(paste("Some of the columns of the",mcall,"data frame have non-numeric variables.\n"))
  }
  
  unvalidobs <- which(apply(MicDtDF,1,function(v) all(!is.finite(v))))
  nunvalid <- length(unvalidobs) 
  if (nunvalid>0) {
    MicDtDF <- MicDtDF[-unvalidobs,]
    agrby <- agrby[-unvalidobs]
    string2 <- paste("rows of the",mcall,"data frame were dropped because they only included non-valid (non finite or missing values) observations.\n")
    if (nunvalid<=10) warning(paste("The",paste(row.names(MicDtDF)[unvalidobs],collapse=" "),string2))
    else warning(paste(nunvalid,string2,collapse=" "))
  }

  if (!is.factor(agrby)) stop("Argument agrby is not a factor\n")
  globaln <- nrow(MicDtDF)
  if (length(agrby)!=globaln) stop("Size of the agrby argument does not agree with the number of rows in the MicDtDF data frame.\n") 
  if ( agrcrt[1]!="minmax" && (!is.numeric(agrcrt) || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
    stop(paste("Wrong value for the agrcrt argument\n( it should be either the string minmax or a two-dim vector",
               "\nof a prob. value for the lower percentile, followed by the prob. value for the upper percentile - \nex:c(0.05,0.95) ).\n")) 
  
  #Trim microdata, if applicable
  if (agrcrt[1] != "minmax") {
    MicDtDF_trim <- MicDtDF
    split_idx <- split(seq_len(nrow(MicDtDF)), agrby)
    for (grp in names(split_idx)) {
      idx <- split_idx[[grp]]
      for (col in seq_len(ncol(MicDtDF))) {
        x <- MicDtDF[idx, col]
        if (!all(is.na(x))) {
          q <- quantile(x, probs = agrcrt, na.rm = TRUE)
          sel <- which(x < q[1] | x > q[2])
          if (length(sel) > 0) MicDtDF_trim[idx[sel], col] <- NA
        }
      }
    }
    MicDtDF <- MicDtDF_trim
  }

  if (length(unique(agrby))!=length(levels(agrby)))  agrby <- factor(agrby)
  grplvls <- levels(agrby)
  nvar <- ncol(MicDtDF)

  # logical vector: TRUE = group is valid; FALSE = group has at least one variable all NA
  keep_group <- sapply(grplvls, function(g) {
    rind <- which(agrby == g)
    all_na_in_any_var <- any(sapply(1:nvar, function(c) all(is.na(MicDtDF[rind, c]))))
    !all_na_in_any_var  # keep if FALSE
  })

  dropped_groups <- grplvls[!keep_group]
  if (length(dropped_groups) > 0) {
    warning(sprintf(
      "Removed %d groups with at least one variable fully NA: %s",
      length(dropped_groups),
      paste(dropped_groups, collapse = ", ")
    ))
    
    # valid group levels
    valid_grplvls <- grplvls[keep_group]

    # subset data and grouping factor
    valid_idx <- agrby %in% valid_grplvls
    MicDtDF <- MicDtDF[valid_idx, , drop = FALSE]
    agrby   <- droplevels(agrby[valid_idx])
    grplvls <- valid_grplvls
    nvar <- ncol(MicDtDF)
  }

  ngrps <- length(grplvls)
  bndsDF <- as.data.frame(matrix(NA_real_, nrow = ngrps, ncol = 2 * nvar))
  NbMicroUnits <- integer(ngrps)
  for (r in 1:ngrps) { 
    grp <- grplvls[r]
    rind <- which(agrby==grp)
    NbMicroUnits[r] <- length(rind)
    for (c in 1:nvar) {
      bndsDF[r,c] <- min(MicDtDF[rind,c], na.rm = TRUE)
      bndsDF[r,nvar+c] <- max(MicDtDF[rind,c], na.rm = TRUE)
    }
  }
  rownames(bndsDF)<-grplvls
  if (!identical(LatentCase, c("U_id_symmetric","U_id","General"))){LatentCase <- match.arg(LatentCase)}
    else {
      if (is.null(LatentParam)||length(LatentParam)==1) {LatentCase <- "U_id_symmetric"}
      else if (assertthat::is.number(LatentParam[[1]])) {LatentCase <- "U_id"}
      else {LatentCase <- "General"}
    }
  if(LatentCase=="General"&&identical(LatentDist, c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"))) LatentDist<-"KDE"
  Umicro<-get_latent_var(MicDtDF,bndsDF,agrby,agrlevels=grplvls,Seq="AllLb_AllUb")
  res <- intData(bndsDF,Seq="AllLb_AllUb",LatentParam,LatentCase,LatentDist,TriangParam,BetaParam.a,BetaParam.b,Umicro,estimate.DistParam,VarNames=names(MicDtDF),ObsNames=grplvls)
  DegInT <- which(apply(res@Ranges,1,function(v) any(!is.finite(v))))
  nDegInT <- length(DegInT)
  if (nDegInT>0) {
    if (nDegInT==res@NObs) {
      warning("No Idata object was created because all units had some degenerate intervals")
      return(NULL)
    }
    if (nDegInT<10) {
      if (nDegInT==1) {
        wmsg <- paste("Data unit",res@ObsNames[DegInT],"was eliminated because it lead to some degenerate intervals")
      } else {
        wmsg <- paste(
          "Data units",paste(res@ObsNames[DegInT],collapse=", "),"were eliminated because they lead to some degenerate intervals",sep="\n"
        )
      }  
    } else {
      wmsg <- paste(nDegInT,"data units were eliminated because they lead to some degenerate intervals")
    }
    warning(wmsg)
    res <- res[-DegInT,]
    res@NbMicroUnits <- NbMicroUnits[-DegInT]
  } else {
    res@NbMicroUnits <- NbMicroUnits
  }  
  names(res@NbMicroUnits) <- res@ObsNames
  res
}
