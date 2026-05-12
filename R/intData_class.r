#' Interval Data Class
#'
#' A class to represent interval data.
#' 
#' @slot Centers A data frame of centers of the intervals.
#' @slot Ranges A data frame of ranges of the intervals.
#' @slot LatentParam A list with the parameters of the latent variables.
#' @slot LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @slot LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"},\code{"Triang"},\code{"TNorm"},\code{"InvTri"},\code{"Beta"},\code{"KDE"},\code{"Degenerated"}), if not, it is a vector with the distribution for each variable.
#' @slot ObsNames A character vector of observation names.
#' @slot VarNames A character vector of variable names.
#' @slot NObs A numeric value indicating the number of observations.
#' @slot NIVar A numeric value indicating the number of interval variables.
#' @slot NbMicroUnits An integer indicating the number of micro units.
#' @import methods
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
#' @references Adapted from package \code{MAINT.Data} (\url{https://cran.r-project.org/package=MAINT.Data}).
#' @export
setClass("intData",slots=c(
  Centers="data.frame",
  Ranges="data.frame",
  LatentParam="list",
  LatentCase="character",
  LatentDist="character",
  ObsNames="character",
  VarNames="character",
  NObs="numeric",
  NIVar="numeric",
  NbMicroUnits="integer"
))

#' Summary Interval Data Class
#'
#' A class to represent the summary of interval data.
#'
#' @slot Centersumar A table summarizing the centers.
#' @slot Rngsumar A table summarizing the ranges.
#' @import methods
#' @export
setClass("summaryintData",slots=c(
  Centersumar="table",
  Rngsumar="table"
))

setGeneric("nrow")
setGeneric("ncol")
setGeneric("rownames")
setGeneric("row.names")
setGeneric("colnames")
setGeneric("var")
setGeneric("plot",signature=c("x","y"))
setGeneric("summary",signature="object")
setGeneric("head",package="utils",signature="x")
setGeneric("tail",package="utils",signature="x")

#' @export
#' @rdname Centers
setGeneric("Centers",function(Sdt) standardGeneric("Centers"))

#' @export
#' @rdname LogRanges
setGeneric("LogRanges",function(Sdt) standardGeneric("LogRanges"))

#' @export
#' @rdname Ranges
setGeneric("Ranges",function(Sdt) standardGeneric("Ranges"))

#' @export
#' @rdname LowerBounds
setGeneric("LowerBounds",function(Sdt) standardGeneric("LowerBounds"))

#' @export
#' @rdname UpperBounds
setGeneric("UpperBounds",function(Sdt) standardGeneric("UpperBounds"))

#' @export
#' @rdname NbMicroUnits
setGeneric("NbMicroUnits",function(x) standardGeneric("NbMicroUnits"))

#' @export
#' @rdname LatentParam
setGeneric("LatentParam",function(Sdt) standardGeneric("LatentParam"))

#' @export
#' @rdname LatentCase
setGeneric("LatentCase",function(Sdt) standardGeneric("LatentCase"))

#' @export
#' @rdname LatentDist
setGeneric("LatentDist",function(Sdt) standardGeneric("LatentDist"))

#' Interval Data Constructor
#'
#' Constructs an interval data object.
#'
#' @param Data A data frame or matrix containing the data.
#' @param Seq Format of macrodata if it is a data frame or matrix. Available options are:
#' \itemize{
#'   \item \code{"AllLb_AllUb"}: All lower bounds followed by all upper bounds, in the same variable order.
#'   \item \code{"AllCen_AllRng"}: All Centers followed by all Ranges, in the same variable order.
#'   \item \code{"LbUb_VarbyVar"}: Lower bounds followed by upper bounds, variable by variable.
#'   \item \code{"CenRng_VarbyVar"}: Centers followed by Ranges, variable by variable.
#' }
#' @param LatentParam A list with the parameters of the latent variables.
#' @param LatentCase A string specifying which of the three scenarios applies to the latent variables:
#' \itemize{
#'   \item \code{"General"}: The case where the latent variables do not have any nice properties.
#'   \item \code{"U_id"}: The case where the latent variables are identically distributed.
#'   \item \code{"U_id_symmetric"}: The case where the latent variables are identically distributed and symmetric.
#' }
#' Defaults to \code{"U_id_symmetric"}.
#' @param LatentDist A string or vector of strings specifying the distribution(s) of the latent variables. If the variables are identically distributed it can be one of (\code{"Unif"},\code{"Triang"},\code{"TNorm"},\code{"InvTri"},\code{"Beta"},\code{"KDE"},\code{"Degenerated"}), if not a vector must be provided with the distribution for each variable.
#' @param TriangParam Mode of the triangular distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{0}.
#' @param BetaParam.a Parameter alpha of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param BetaParam.b Parameter beta of the Beta distribution. If the latent variables are identically distributed, it is only necessary to provide a number, if not a vector is needed.
#' The default is \code{1}.
#' @param Umicro Latent microdata observations. Needed if \code{LatentDist="KDE"} or \code{estimate.DistParam=TRUE}.
#' @param estimate.DistParam Logical parameter indicating if estimation of the parameters of the latent distributions should be performed. Can only be set to TRUE if \code{LatentCase="General"}.
#' The default is \code{FALSE}.
#' @param VarNames A character vector of variable names.
#' @param ObsNames A character vector of observation names.
#' @param NbMicroUnits An integer specifying the number of micro units.
#' 
#' @importFrom assertthat is.number
#' 
#' @return An object of class \linkS4class{intData}.
#' 
#' @references Oliveira, M. R., Pinheiro, D., & Oliveira, L. (2025). 
#' Location and association measures for interval-valued data based on Mallows' distance. 
#' arXiv preprint arXiv:2407.05105. \url{https://arxiv.org/abs/2407.05105}
#' @references Adapted from package \code{MAINT.Data} (\url{https://cran.r-project.org/package=MAINT.Data}).
#' 
#' @export
intData <- function(Data,
                    Seq=c("AllLb_AllUb","AllCen_AllRng","LbUb_VarbyVar","CenRng_VarbyVar"),
                    LatentParam=NULL,
                    LatentCase=c("U_id_symmetric","U_id","General"),
                    LatentDist=c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated"),
                    TriangParam=0,
                    BetaParam.a=1,
                    BetaParam.b=1,
                    Umicro=NULL,
                    estimate.DistParam=FALSE,
                    VarNames=NULL,
                    ObsNames=row.names(Data),
                    NbMicroUnits=integer(0)){

  if ( !is.data.frame(Data) && !is.matrix(Data) ) stop("First argument of intData must be a data frame or a matrix\n")
  if (!is.integer(NbMicroUnits)) {
    unitnames <- names(NbMicroUnits)
    NbMicroUnits <- as.integer(NbMicroUnits)
    names(NbMicroUnits) <- unitnames
  }  
  
  p <- ncol(Data)  # Total number of Interval variable bounds
  q <- p/2	 # Number of Interval variables
  if (floor(q) != q) stop("Number of columns of Data ( =",p,") must be an even number\n")
  Seq <- match.arg(Seq)
  if (  (Seq == "LbUb_VarbyVar") || (Seq == "AllLb_AllUb") )
  {
    if (Seq == "LbUb_VarbyVar") { Lbnd <- Data[,2*(0:(q-1))+1] ; Ubnd <- Data[,2*(1:q)] }
    if (Seq == "AllLb_AllUb")   { Lbnd <- Data[,1:q] ; Ubnd <- Data[,(q+1):p] }
    Centers <- (Lbnd+Ubnd)/2
    Ranges <- (Ubnd-Lbnd)
    if (any(is.na(Ranges))) stop("Invalid data")
  } else {
    if (Seq == "CenRng_VarbyVar") { Centers <- Data[,2*(0:(q-1))+1] ; Ranges <- Data[,2*(1:q)] }
    if (Seq == "AllCen_AllRng")   { Centers <- Data[,1:q] ; Ranges <- Data[,(q+1):p] }
  }
  if (is.null(VarNames)) VarNames <- paste("I",1:q,sep="")
  if (!is.data.frame(Centers)) Centers <- as.data.frame(Centers)
  if (!is.data.frame(Ranges)) Ranges <- as.data.frame(Ranges)
  names(Centers) <- paste(VarNames,".Centers",sep="")
  names(Ranges) <- paste(VarNames,".Ranges",sep="")
  if (is.null(ObsNames)) ObsNames <- as.character(seq_len(nrow(Centers)))
  rownames(Centers) <- rownames(Ranges) <- ObsNames

  if (!identical(LatentCase, c("U_id_symmetric","U_id","General"))){LatentCase <- match.arg(LatentCase)}
    else {
      if (is.null(LatentParam)||length(LatentParam)==1) {LatentCase <- "U_id_symmetric"}
      else if (assertthat::is.number(LatentParam[[1]])) {LatentCase <- "U_id"}
      else {LatentCase <- "General"}
    }
  
  if(LatentCase!="General"&&length(unique(TriangParam)) > 1) stop("Error: For different TriangParam for each variable, LatentCase must be 'General'.")
  if(LatentCase!="General"&&length(unique(BetaParam.a)) > 1) stop("Error: For different BetaParam.a for each variable, LatentCase must be 'General'.")

  if (is.null(LatentParam)) LatentParam<-get_latent_param(LatentCase,LatentDist,TriangParam,BetaParam.a,BetaParam.b,Umicro,q,estimate.DistParam)[[1]]

  new("intData",Centers=Centers,Ranges=Ranges,LatentParam=LatentParam,LatentCase=LatentCase,LatentDist=LatentDist,
      ObsNames=ObsNames,VarNames=VarNames,NObs=nrow(Centers),NIVar=q,NbMicroUnits=NbMicroUnits)
}

#' Summary Method for \linkS4class{intData}
#'
#' @param object An object of class \linkS4class{intData}.
#' @return An object of class \code{summaryintData}.
#' @import methods
#' @export
#' @rdname summary
#' @aliases summary,intData-method
setMethod("summary",
  signature(object = "intData"),
  function (object) 
  {
    show(new("summaryintData",Centersumar=summary(object@Centers),Rngsumar=summary(object@Ranges)))
  }
)

#' Show Method for \linkS4class{intData}
#'
#' @param object An object of class \linkS4class{intData}.
#' @return The object itself, returned invisibly. Called for its side effects (printing).
#' @import methods
#' @export
#' @rdname show
#' @aliases show,intData-method
setMethod("show",
  signature(object = "intData"),
  function (object) 
  {
    printrow <- function(Bnds,NIVar) 
    { 
      cat(Bnds[1],"  ")
      for (j in 2:(NIVar+1))
        cat("[",format(Bnds[j],width=8,digits=5,justify="centre"),", ",
          format(Bnds[NIVar+j],width=8,digits=5,justify="centre"),"]  ",sep="") 
      cat("\n")
    }

    HalfRange <- object@Ranges/2
    LB <- object@Centers - HalfRange
    UB <- object@Centers + HalfRange
    lobsname <- max(nchar(object@ObsNames))
    flength <- max(nchar(object@VarNames),nchar(format(LB[1,1],width=8,digits=5))+nchar(format(UB[1,1],width=8,digits=5))) + 6 
    for (j in 1:object@NIVar) {
      if (j>1) {
        nspaces <- flength-nchar(object@VarNames[j])
      } else {
        nspaces <- lobsname+ceiling(flength/2)-ceiling(nchar(object@VarNames[1])/2)+3
      }  
      cat(rep(" ",nspaces),object@VarNames[j],sep="" )
    }  
    cat("\n") 
    apply(cbind(format(object@ObsNames,width=lobsname),LB,UB),1,printrow,NIVar=object@NIVar)
    invisible(object)
  }
)

#' Number of Rows Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return The number of rows.
#' @import methods
#' @export
#' @rdname nrow
#' @aliases nrow,intData-method
setMethod("nrow",signature(x = "intData"),function(x) x@NObs)

#' Number of Columns Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return The number of columns.
#' @import methods
#' @export
#' @rdname ncol
#' @aliases ncol,intData-method
setMethod("ncol",signature(x = "intData"),function(x) x@NIVar)

#' Dimensions Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return A vector of the number of rows and columns.
#' @import methods
#' @export
#' @rdname dim
#' @aliases dim,intData-method
setMethod("dim",signature(x = "intData"),function(x) c(nrow(x),ncol(x)))

#' Row Names Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return A character vector of row names.
#' @import methods
#' @export
#' @rdname rownames
#' @aliases rownames,intData-method
setMethod("rownames",signature(x = "intData"),function(x) x@ObsNames)

#' Row.Names Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return A character vector of row names.
#' @import methods
#' @export
#' @rdname row.names
#' @aliases row.names,intData-method
setMethod("row.names",signature(x = "intData"),function(x) x@ObsNames)

#' Column Names Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return A character vector of column names.
#' @import methods
#' @export
#' @rdname colnames
#' @aliases colnames,intData-method
setMethod("colnames",signature(x = "intData"),function(x) x@VarNames)

#' Variable Names Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return A character vector of variable names.
#' @import methods
#' @export
#' @rdname names
#' @aliases names,intData-method
setMethod("names",signature(x = "intData"),function(x) x@VarNames)

#' Centers Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A \code{data.frame} containing the centers of the intervals.
#' @import methods
#' @export
#' @rdname Centers
setMethod("Centers",signature(Sdt = "intData"),function(Sdt) Sdt@Centers)

#' LogRanges Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A \code{data.frame} containing the logarithms of the ranges.
#' @import methods
#' @export
#' @rdname LogRanges
setMethod("LogRanges",signature(Sdt = "intData"),function(Sdt) log(Sdt@Ranges))

#' Ranges Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A \code{data.frame} containing the ranges of the intervals.
#' @import methods
#' @export
#' @rdname Ranges
setMethod("Ranges",signature(Sdt = "intData"),function(Sdt) Sdt@Ranges)

#' Lower Bounds Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A \code{data.frame} containing the lower bounds of the intervals.
#' @import methods
#' @export
#' @rdname LowerBounds
setMethod("LowerBounds",signature(Sdt = "intData"),function(Sdt){
  bounds <- Sdt@Centers-0.5*Sdt@Ranges
  names(bounds) <- paste(Sdt@VarNames,".Lbnd",sep="")
  return(bounds)
})

#' Upper Bounds Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A \code{data.frame} containing the upper bounds of the intervals.
#' @import methods
#' @export
#' @rdname UpperBounds
setMethod("UpperBounds",signature(Sdt = "intData"),function(Sdt){
  bounds <- Sdt@Centers+0.5*Sdt@Ranges
  names(bounds) <- paste(Sdt@VarNames,".Ubnd",sep="")
  return(bounds)
})

#' Latent Parameters Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A list with the latent parameters.
#' @import methods
#' @export
#' @rdname LatentParam
setMethod("LatentParam",signature(Sdt = "intData"),function(Sdt) Sdt@LatentParam)

#' Latent Case Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A character with the latent case.
#' @import methods
#' @export
#' @rdname LatentCase
setMethod("LatentCase",signature(Sdt = "intData"),function(Sdt) Sdt@LatentCase)

#' Latent Distribution Method for \linkS4class{intData}
#'
#' @param Sdt An object of class \linkS4class{intData}.
#' @return A character with the latent distribution(s).
#' @import methods
#' @export
#' @rdname LatentDist
setMethod("LatentDist",signature(Sdt = "intData"),function(Sdt) Sdt@LatentDist)

#' Number of Micro Units Method for \linkS4class{intData}
#'
#' @param x An object of class \linkS4class{intData}.
#' @return An integer specifying the number of micro units.
#' @import methods
#' @export
#' @rdname NbMicroUnits
setMethod("NbMicroUnits",
  signature(x = "intData"),
  function(x) {
    if (length(x@NbMicroUnits)==0) return(NULL)
    x@NbMicroUnits
  }
)

#' Head Method for \linkS4class{intData}
#'
#' Returns the first \code{n} rows of an \linkS4class{intData} object.
#'
#' @param x An \linkS4class{intData} object.
#' @param n The number of rows to return.
#' @return A subset of the \linkS4class{intData} object.
#' @import methods
#' @export
#' @rdname head
#' @aliases head,intData-method
setMethod("head",
  signature(x = "intData"),
  function (x,n=min(nrow(x),6L)) 
  {
    
    if (n>0) {
      x[1:n,1:x@NIVar] 
    } else {
      x[-1:n,1:x@NIVar]
    }
  }
)

#' Tail Method for \linkS4class{intData}
#'
#' Returns the last \code{n} rows of an \linkS4class{intData} object.
#'
#' @param x An \linkS4class{intData} object.
#' @param n The number of rows to return.
#' @return A subset of the \linkS4class{intData} object.
#' @import methods
#' @export
#' @rdname tail
#' @aliases tail,intData-method
setMethod("tail",
  signature(x = "intData"),
  function (x,n=min(nrow(x),6L)) 
  {
    if (n>0) {
      x[(x@NObs-n+1):x@NObs,1:x@NIVar] 
    } else {
      x[(-x@NObs-n-1):-x@NObs,1:x@NIVar]
    }
  }
)

#' Plot Method for Two \linkS4class{intData} Objects
#'
#' Plots one \linkS4class{intData} object against another, with options to visualize the intervals as crosses or rectangles.
#'
#' @param x An \linkS4class{intData} object to plot on the x-axis.
#' @param y An \linkS4class{intData} object to plot on the y-axis.
#' @param type The type of plot to generate: "crosses" or "rectangles" or "crosses2". Default is "crosses".
#' @param append Logical, if \code{TRUE}, the plot is added to the current plot.
#' @param palette A vector with colors for each observation.
#' @param ... Additional graphical parameters.
#' @return A plot showing the relationship between the two \linkS4class{intData} objects.
#' @import methods
#' @importFrom graphics plot.default lines rect
#' @export
#' @rdname plot
#' @aliases plot,intData,intData-method
setMethod("plot",
  signature(x = "intData",y = "intData"),
  function(x, y, type=c("crosses","rectangles","crosses2"), append=FALSE, palette=rainbow(x@NObs), ...)
  {
    if (x@NIVar > 1) stop("Currently intData method plot can plot only one integer variable on the horizontal axis\n")
    if (y@NIVar > 1) stop("Currently intData method plot can plot only one integer variable on the vertical axis\n")
    if (x@NObs != y@NObs) stop("Arguments x and y have a different number of elements\n")

    type <- match.arg(type)

    dotarguments <- match.call(expand.dots=FALSE)$...

    if (is.null(dotarguments$main)) {
      dotarguments$main <- paste(y@VarNames,"vs.",x@VarNames)
    }
    if (is.null(dotarguments$xlab)) {
      dotarguments$xlab <- x@VarNames
    }
    if (is.null(dotarguments$ylab)) {
      dotarguments$ylab <- y@VarNames
    }


    HlfRngx <- x@Ranges[,1]/2
    Lbx <- x@Centers[,1] - HlfRngx
    Ubx <- x@Centers[,1] + HlfRngx
    HlfRngy <- y@Ranges[,1]/2
    Lby <- y@Centers[,1] - HlfRngy
    Uby <- y@Centers[,1] + HlfRngy

    if (is.null(dotarguments$xlim)) {
      minvx <- min(Lbx)
      maxvx <-max(Ubx)
      dotarguments$xlim <- c(0.95*minvx,1.05*maxvx)
    }
    if (is.null(dotarguments$ylim)) {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$ylim <- c(0.95*minvy,1.05*maxvy)
    }

    if (!append) {
      do.call("plot.default",c(list(x=0.,y=0.,type="n"),dotarguments))
    }

    if (type=="crosses") {
      for(i in 1:x@NObs){
        lines(c(Lbx[i],Ubx[i]),c((Lby[i]+Uby[i])/2,(Lby[i]+Uby[i])/2),col=palette[i],lwd=1,...)
        lines(c((Lbx[i]+Ubx[i])/2,(Lbx[i]+Ubx[i])/2),c(Lby[i],Uby[i]),col=palette[i],lwd=1,...)
      }
    } else if (type=="rectangles") {
      for (i in 1:x@NObs) rect(Lbx[i],Lby[i],Ubx[i],Uby[i],lty=1,border=palette[i], ...)
    }else if (type=="crosses2"){
      for(i in 1:x@NObs){
        lines(c(Lbx[i],Ubx[i]),c(Lby[i],Uby[i]),col=palette[i],lwd=1,...)
        lines(c(Lbx[i],Ubx[i]),c(Uby[i],Lby[i]),col=palette[i],lwd=1,...)
      }
    }
  }
)

#' Plot Method for a Single \linkS4class{intData} Object
#'
#' Plots a single \linkS4class{intData} object, either in a vertical or horizontal layout.
#'
#' @param x An \linkS4class{intData} object.
#' @param casen A vector specifying the case numbers to plot. Default is \code{NULL}.
#' @param layout The layout of the plot: "vertical" or "horizontal".
#' @param append Logical, if \code{TRUE}, the plot is added to the current plot.
#' @param ... Additional graphical parameters.
#' @return A plot showing the intervals of the \linkS4class{intData} object.
#' @import methods
#' @importFrom stats runif
#' @importFrom graphics plot.default segments
#' @export
#' @rdname plot
#' @aliases plot,intData,missing-method
setMethod("plot",
  signature(x = "intData",y = "missing"),
  function(x, casen=NULL, layout=c("vertical","horizontal"), append=FALSE,  ...)
  {
    if (x@NIVar > 1) {
      if (x@NIVar==2) {
        plot.default(x[,1],x[,2],...)
        return()
      } else {
        stop("Currently method plot can only plot at most two interval variables simultaneously\n")
      }
    }

    layout <- match.arg(layout)

    dotarguments <- match.call(expand.dots=FALSE)$...
    if (is.null(dotarguments$main)) {
      dotarguments$main <- x@VarNames
    }
    if (is.null(dotarguments$ylab)) {
      if (layout=="vertical") dotarguments$ylab <- x@VarNames
      else dotarguments$ylab <- "Case numbers"
    }
    if (is.null(dotarguments$xlab)) {
        if (layout=="vertical") dotarguments$xlab <- "Case numbers"
        else dotarguments$xlab <- x@VarNames
    }

    HlfRngy <- x@Ranges[,1]/2
    Lby <- x@Centers[,1] - HlfRngy
    Uby <- x@Centers[,1] + HlfRngy
    xcord <- 1:x@NObs
    if (!append) xcord <- 1:x@NObs
    else xcord <- 1:x@NObs + runif(x@NObs,-0.5,0.5)

    if (is.null(dotarguments$ylim) && layout=="vertical") {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$ylim <- c(0.95*minvy,1.05*maxvy)
    }
    if (is.null(dotarguments$xlim) && layout=="horizontal") {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$xlim <- c(0.95*minvy,1.05*maxvy)
    }

    if (is.null(casen)) {
      casen <- c(0,x@NObs+1)
    } else {
      casen <- factor(eval(casen))
    }
    if (layout=="vertical") {
      if (!append) {
        do.call("plot.default",c(list(x=casen,y=rep(0.,length(casen)),type="n"),dotarguments))
      }
      do.call("segments",c(list(x0=xcord,y0=Lby,x1=xcord,y1=Uby),dotarguments))
    } else  {
      if (!append) {
        do.call("plot.default",c(list(x=rep(0.,length(casen)),y=casen,type="n"),dotarguments))
      }
      do.call("segments",c(list(y0=xcord,x0=Lby,y1=xcord,x1=Uby),dotarguments))
    }
  }
)

#' Show Method for Summary \linkS4class{intData}
#'
#' @param object An object of class \code{summaryintData}.
#' @import methods
#' @export
#' @rdname show
#' @aliases show,summaryintData-method
setMethod("show",
  signature(object = "summaryintData"),
  function(object)
  {
    cat("Centers summary:\n") ; print(object@Centersumar)
    cat("Ranges summary:\n") ; print(object@Rngsumar)
    invisible(object)
  }
)

#' Print Method for Summary \linkS4class{intData}
#'
#' @param x An object of class \code{summaryintData}.
#' @param ... Additional arguments passed to print.
#' @return The object itself, returned invisibly. Called for its side effects (printing).
#' @import methods
#' @export
#' @rdname print
#' @aliases print,summaryintData-method
setMethod("print", signature(x="summaryintData"), function(x,...) invisible(x) )

#-----Standard operators for intData objects-----

#---Indexing and assignement---

#' Subset an \linkS4class{intData} Object
#'
#' Extract a subset of rows and columns from an \linkS4class{intData} object.
#'
#' @param x An \linkS4class{intData} object.
#' @param i Row indices or names to subset. Defaults to all rows.
#' @param j Column indices or names to subset. Defaults to all columns.
#' @param ... Additional arguments (not used).
#' @param drop Logical, passed to the underlying \code{[}. Defaults to \code{TRUE}.
#' @return An \linkS4class{intData} object containing the specified subset of rows and columns.
#' @import methods
#' @export
#' @rdname S4-Extract-methods
#' @aliases [,intData,ANY,ANY-method
setMethod("[",
  signature(x='intData',i='ANY',j='ANY'),
  function(x,i,j,...,drop=TRUE)
  {
    if (missing(i)) i <- seq_len(nrow(x))
    if (missing(j)) j <- seq_len(ncol(x))
    if (is.character(i)) i <- match(i, x@ObsNames)
    if (is.character(j)) j <- match(j, x@VarNames)
    latent_dist <- ifelse(length(x@LatentDist)>1, x@LatentDist[j], x@LatentDist)
    if (x@LatentCase == "General") {
      latent_param <- list(
        x@LatentParam[[1]][j, j, drop = FALSE],
        x@LatentParam[[2]][j, j, drop = FALSE]
      )
    } else {
      latent_param <- x@LatentParam
    }
    intData(cbind(x@Centers[i,j,drop=FALSE],x@Ranges[i,j,drop=FALSE]),
      Seq="AllCen_AllRng",VarNames=x@VarNames[j],ObsNames=x@ObsNames[i],
      LatentParam=latent_param,LatentCase=x@LatentCase,LatentDist=latent_dist)
  }
)

#---Comparison---

#' Equality Comparison for \linkS4class{intData} Objects
#'
#' Compare two \linkS4class{intData} objects for equality.
#'
#' @param e1 First \linkS4class{intData} object.
#' @param e2 Second \linkS4class{intData} object.
#' @return A logical matrix indicating which elements are equal between the two \linkS4class{intData} objects.
#' @import methods
#' @export
#' @rdname comparison-methods
#' @aliases ==
#' @aliases ==.intData
#' @aliases ==,intData,intData-method
setMethod("==",
  signature(e1 = "intData",e2 = "intData"),
  function(e1,e2)
  {
    CompIvalue <- function(Ival) Ival[1,1] == Ival[1,2] && Ival[2,1] == Ival[2,2]

    if (!is(e2,"intData")) 
      stop("Trying to compare an intData object with an object of a diferent type\n")
    if ( e1@NObs != e2@NObs || e1@NIVar != e2@NIVar )
      stop("== only defined for equally-sized intData objects\n")
    if (e1@LatentCase != e2@LatentCase) return(FALSE)
    if (any(e1@LatentDist!= e2@LatentDist)) return(FALSE)
    if (!identical(e1@LatentParam,e2@LatentParam)) return(FALSE)
    TmpArray <- array(dim=c(e1@NObs,e1@NIVar,2,2))
    for (j in 1:e1@NIVar)  {
      TmpArray[,j,,1] <- cbind(e1@Centers[,j],e1@Ranges[,j])
      TmpArray[,j,,2] <- cbind(e2@Centers[,j],e2@Ranges[,j])
    }

    apply(TmpArray,c(1,2),CompIvalue)
  }
)

#' Inequality Comparison for \linkS4class{intData}
#'
#' Compare two \linkS4class{intData} objects for inequality.
#' 
#' @param e1 An \linkS4class{intData} object.
#' @param e2 An \linkS4class{intData} object.
#' @return A logical matrix indicating element-wise inequality of the two \linkS4class{intData} objects.
#' @import methods
#' @export
#' @rdname comparison-methods
#' @aliases !=
#' @aliases !=.intData
#' @aliases !=,intData,intData-method
setMethod("!=",
  signature(e1 = "intData",e2 = "intData"),
  function(e1,e2)  
  {
    if (!is(e2,"intData"))
      stop("Trying to compare an intData object with an object of a diferent type\n")
    if ( e1@NObs != e2@NObs || e1@NIVar != e2@NIVar )
      stop("!= only defined for equally-sized intData objects\n")
    !(e1==e2)
  }
)
