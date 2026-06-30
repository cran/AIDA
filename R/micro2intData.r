#' Aggregate Microdata into Interval-Valued Data
#'
#' Aggregates microdata from a data frame into interval-valued data using various criteria and latent distribution settings.
#'
#' @param microdata A data frame containing the microdata. All columns should be numeric.
#' @param agrby A factor used to specify the grouping of the microdata for aggregation.
#' @param agrcrt A string or numeric vector of length 2 specifying the aggregation criterion. The default is \code{"minmax"}, which takes the minimum and maximum values for each variable. If a numeric vector is provided, it should specify the lower and upper percentiles for aggregation (e.g., \code{c(0.05, 0.95)}).
#' @param LatentParam (Optional) A list with the parameters of the latent variables. 
#'  Expects a list with a single number if `LatentCase` is `"U_id_symmetric"`, a list of two numbers if `LatentCase` is `"U_id"`, and a list of two matrices if `LatentCase` is `"General"`.
#' @param LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"}, \code{"Triang"}, \code{"TNorm"}, \code{"InvTri"}, \code{"Beta"}, \code{"KDE"}, \code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' The default is \code{"Unif"} if \code{LatentCase="U_id_symmetric"}, and \code{"KDE"} if \code{LatentCase="General"}.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param estimate.DistParam Logical parameter indicating if estimation of the parameters of the latent distributions should be performed. Can only be set to TRUE if \code{LatentCase="General"}.
#' The default is \code{FALSE}.
#'
#' @return An \code{\linkS4class{intData}} object containing the aggregated interval-valued data, or \code{NULL} if all units lead to degenerate intervals.
#'
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
#' 
#' # Define grouping variable for microdata aggregation
#' credit_agrby <- factor(paste(CreditCard_microdata$Name, CreditCard_microdata$Month, sep = "_"))
#' 
#' # Create intData object by aggregating microdata using the default minmax criterion 
#' # and using KDE for estimation of the latent distribution in the general case
#' credit_agr <- micro2intData(CreditCard_microdata[,3:7],
#'                             agrby = credit_agrby,
#'                             LatentCase = "General")
#' 
#' @export
micro2intData <- function(microdata,
                          agrby,
                          agrcrt="minmax",
                          LatentParam=NULL,
                          LatentCase=c("U_id_symmetric","U_id","General"),
                          LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                          TriangParam=0,
                          BetaParam.a=1,
                          BetaParam.b=1,
                          estimate.DistParam=FALSE){
  mcall <- match.call()$microdata
  if (length(mcall) > 1) mcall <- "microdata"
  if (!(is.data.frame(microdata))) stop("First argument of AgMicroData must be a data frame.\n")
  if (!is.data.frame(microdata)) microdata <- as.data.frame(microdata)
  if (any(!sapply(seq_len(ncol(microdata)),function(ind) is.numeric(microdata[,ind])))){  
    stop(paste("Some of the columns of the",mcall,"data frame have non-numeric variables.\n"))
  }
  
  unvalidobs <- which(apply(microdata,1,function(v) all(!is.finite(v))))
  nunvalid <- length(unvalidobs) 
  if (nunvalid>0) {
    microdata <- microdata[-unvalidobs,]
    agrby <- agrby[-unvalidobs]
    string2 <- paste("rows of the",mcall,"data frame were dropped because they only included non-valid (non finite or missing values) observations.\n")
    if (nunvalid<=10) warning(paste("The",paste(row.names(microdata)[unvalidobs],collapse=" "),string2))
    else warning(paste(nunvalid,string2,collapse=" "))
  }

  if (!is.factor(agrby)) stop("Argument agrby is not a factor\n")
  globaln <- nrow(microdata)
  if (length(agrby)!=globaln) stop("Size of the agrby argument does not agree with the number of rows in the microdata data frame.\n") 
  if ( agrcrt[1]!="minmax" && (!is.numeric(agrcrt) || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
    stop(paste("Wrong value for the agrcrt argument\n( it should be either the string minmax or a two-dim vector",
               "\nof a prob. value for the lower percentile, followed by the prob. value for the upper percentile - \nex:c(0.05,0.95) ).\n")) 
  
  # Trim microdata, if applicable
  if (agrcrt[1] != "minmax") {
    microdata_trim <- microdata
    split_idx <- split(seq_len(nrow(microdata)), agrby)
    for (grp in names(split_idx)) {
      idx <- split_idx[[grp]]
      for (col in seq_len(ncol(microdata))) {
        x <- microdata[idx, col]
        if (!all(is.na(x))) {
          q <- quantile(x, probs = agrcrt, na.rm = TRUE)
          sel <- which(x < q[1] | x > q[2])
          if (length(sel) > 0) microdata_trim[idx[sel], col] <- NA
        }
      }
    }
    microdata <- microdata_trim
  }

  if (length(unique(agrby))!=length(levels(agrby)))  agrby <- factor(agrby)
  grplvls <- levels(agrby)
  NIVar <- ncol(microdata)

  # logical vector: TRUE = group is valid; FALSE = group has at least one variable all NA
  keep_group <- sapply(grplvls, function(g) {
    rind <- which(agrby == g)
    all_na_in_any_var <- any(sapply(1:NIVar, function(c) all(is.na(microdata[rind, c]))))
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
    microdata <- microdata[valid_idx, , drop = FALSE]
    agrby   <- droplevels(agrby[valid_idx])
    grplvls <- valid_grplvls
    NIVar <- ncol(microdata)
  }

  ngrps <- length(grplvls)
  bndsDF <- as.data.frame(matrix(NA_real_, nrow = ngrps, ncol = 2 * NIVar))
  NbMicroUnits <- integer(ngrps)
  for (r in 1:ngrps) { 
    grp <- grplvls[r]
    rind <- which(agrby==grp)
    NbMicroUnits[r] <- length(rind)
    for (c in 1:NIVar) {
      bndsDF[r,c] <- min(microdata[rind,c], na.rm = TRUE)
      bndsDF[r,NIVar+c] <- max(microdata[rind,c], na.rm = TRUE)
    }
  }
  rownames(bndsDF) <- grplvls
  
  Umicro <- get_latent_var(microdata,bndsDF,agrby,agrlevels=grplvls,Seq="AllLb_AllUb")
  res <- intData(bndsDF,Seq="AllLb_AllUb",LatentParam,LatentCase,LatentDist,TriangParam,BetaParam.a,BetaParam.b,Umicro,estimate.DistParam,VarNames=names(microdata),ObsNames=grplvls)
  DegInT <- which(apply(res@Ranges,1,function(v) any(v==0)))
  nDegInT <- length(DegInT)
  if (nDegInT>0) {
    if (nDegInT==res@NObs) {
      warning("No intData object was created because all units had some degenerate intervals")
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
