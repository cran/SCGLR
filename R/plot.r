# remove warnings due to ggplot2 syntax strangeness
utils::globalVariables(c("comp","y","label","angle","hjust"))

#' SCGLR generic plot
#' @export
#' @method plot SCGLR
#' @S3method plot SCGLR
#' @param x an object from SCGLR class.
#' @param \dots optional arguments.
#' @param style describes which plot will be drawn. Style "simple" : correlation plot, "covariates" (to add covariates arrows), "circle" (to add
#' the unit circle), "observations" (to add observations), "threshold" (to select covariates whose sum of square correlations with the two components 
#' of the plane exceeds the given threshold), "factor" (to add one the factors),
#' "predictors" (to add linear predictor arrows). Can be given as a vector or list of strings (eg: c("simple","circle"))
#' or comma separated string (eg: "simple,circle"). Style elements can also be abbreviated (eg: "simp,cir").
#' @param threshold a numeric value. All covariates whose sum of square correlations with the two components of the plane lower than this threshold will be ignored.
#' @param plane a size-2 vector (or comma separated string) indicating which components are plotted (eg: c(1,2) or "1,2").
#' @param factor factor to show (default to first one).
#' @param predictors a vector of character to select the linear predictors displayed.
#' @param covariates a vector of character to select the covariates displayed.
#' @param label.offset offset by which labels should be moved from tip of arrows.
#' @param label.auto whether or not the labels should be rotated according to vector angle.
#' @param label.size relative size for labels (fine tuning).
#' @param expand expand factor for windows size, for example to make room for clipped labels.
#' @return an object of class \code{ggplot}.
#' @examples \dontrun{
#' library(SCGLR)
#' 
#' # load sample data
#' data(genus)
#' 
#' # get variable names from dataset
#' n <- names(genus)
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx <- n[-grep("^gen",n)]   # X <- remaining names
#' 
#' # remove "geology" and "surface" from nx
#' # as surface is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#' 
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#' 
#' # define family
#' fam <- rep("poisson",length(ny))
#' 
#' genus.scglr <- scglr(formula=form,data = genus,family=fam, K=4,
#'  offset=genus$surface)
#' 
#' summary(genus.scglr)
#' 
#' barplot(genus.scglr)
#' 
#' plot(genus.scglr)
#' 
#' plot(genus.scglr,style="circle,cov,predictors,fact")
#' 
#' pairs(genus.scglr)
#'
#' } 
plot.SCGLR <- function(x, ..., style="simple", threshold=0.8, plane=c(1,2), factor=NULL, predictors=NULL, covariates=NULL, label.offset=0.01, label.auto=T,label.size=1, expand=1.0) {
  data <- x
  if(class(data)!="SCGLR") 
    stop("This plot function need an SCGLR result")
  
  if(dim(data$compr)[2]<2)
    stop("At least two axes are needed for this kind of plot!")
  
  # process style
  if(is.character(style)) {
    style <- str_trim(unlist(strsplit(style,",")))
  }

  style<-match.arg(style,c("simple","covariates","circle","threshold","observations","factor","predictors"),
                   several.ok=TRUE)
  
  # predefined styles
  if("simple" %in% style)
    style <- c(style,"covariates","circle")
  
  if("threshold" %in% style)
    style <- c(style,"circle")
  
  if(!is.null(factor)) {
    style <- c(style,"factor")
  }

  if(!is.null(predictors)) {
    style <- c(style,"predictors")
  }

  if(!is.null(covariates)) {
    style <- c(style,"covariates")
  }
  
  # clean
  style <- unique(style)
  
  # customize
  dots <- list(...)
  cust <- function(key, def=NULL) {
    if(is.null(dots[[key]])) {
      return(def) 
    } else {
      return(dots[[key]])
    }
  }
  
  # process plan
  if(is.character(plane)) {
    plane <- as.integer(str_trim(unlist(strsplit(plane,","))))
  }
  
  # sanity checking
  if(length(plane) !=2 ) {
    stop("Plane should have two components!")
  }
  if((min(plane)<1) || (max(plane)>ncol(data$compr))) {
    stop("Invalid components for plane!")
  }
  
  # plan
  axis_names <- colnames(data$compr)[plane]
  
  # check factor
  if("factor" %in% style) {
    if(is.null(factor)) {
      if(is.null(data$xFactors))
        stop("No factor in data!")
      factor <- names(data$xFactors)[1]
      warning("No factor given, assuming first one! (",factor,")!")
    } else {
      if(!factor %in% names(data$xFactors))
        stop("Invalid factor!")
    }
  }
  
  # inertia
  inertia <- data$inertia[plane]
  
  # build base plot
  p <- qplot((-1:1)*expand, (-1:1)*expand, geom="blank")+
    coord_fixed()+
    # thicker x unit arrow
    xlab(paste(axis_names[1],"(",round(100*inertia[1],2),"%)")) + 
    geom_hline(yintercept=0)+
    geom_segment(aes(x=-1.1,xend=1.1,y=0,yend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))+
    # thicker y unit arrow
    ylab(paste(axis_names[2],"(",round(100*inertia[2],2),"%)")) + 
    geom_vline(xintercept=0)+
    geom_segment(aes(y=-1.1,yend=1.1,x=0,xend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))
  
  # plot title
  if("observations" %in% style) {
    if("circle" %in% style) {
      p <- p + ggtitle("Mixed individual and \ncorrelation plot\n")
    } else {
      p <- p + ggtitle("Individual plot\n")
    }
  } else {
    if("circle" %in% style) {
      p <- p + ggtitle("Correlation plot\n")      
    }
  }

  # add unit circle
  if(("circle" %in% style) || ("threshold" %in% style))
    p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(fill=NA)),-1,1,-1,1)

  # add threshold circle
  if("threshold" %in% style)
    p <- p + annotation_custom(circleGrob(r=0.5*threshold,gp=gpar(fill=NA,lty="dashed")),-1,1,-1,1)
  
  # add observations
  if("observations" %in% style) {
    obs <- as.data.frame(data$compr[,plane])
    names(obs) <- c("x","y")
    if(("factor" %in% style)&&(cust("observations.factor",FALSE))) {
      obs<-cbind(obs,data$xFactors[factor])
      p <- p + geom_point(
        aes_string(x="x",y="y",color=factor),
        data=obs,
        size=cust("observations.size",1)
      )
    } else {
      p <- p + geom_point(
        aes(x=x,y=y),
        data=obs,
        size=cust("observations.size",1),
        color=cust("observations.color","black")
      )
    }
  }
  
  # add linear predictor arrows
  if("predictors" %in% style) {
    co <- as.data.frame(cor(data$lin.pred, data$compr[,plane]))
    if(!is.null(predictors)) {
      co <- co[rownames(co)%in%predictors,]
    }
    names(co) <- c("x", "y")
    co$norm <- sqrt(co$x^2+co$y^2)
    co$label <- rownames(co)
    
    # adjust label position
    co$angle <- atan2(co$y,co$x)*180/pi
    co$hjust <- ifelse(abs(co$angle)>90,1,0)
    if(label.auto) {
      co$angle <- ifelse(abs(co$angle)>90,co$angle+180,co$angle)
    } else {
      co$angle <- 0
    }
    
    # filter according to given threshold
    if("threshold" %in% style)
      co <- co[co$norm>threshold,]
    
    if(nrow(co)==0) {
      warning("No correlation higher than threshold value ",threshold,"!")
    } else {
      clr <- cust("predictors.color","red")
      if(cust("predictors.arrows",TRUE)) {
        p <- p + geom_segment(
          aes(x=0,y=0,xend=x,yend=y),
          data=co,
          color=cust("predictors.arrows.color",clr),
          arrow=arrow(length=unit(0.02,"npc")))
      } else {
        # reset label position
        co$hjust<-0.5
        co$angle<-0
      }
      if(cust("predictors.labels",TRUE)) {
        p <- p + geom_text(
          aes(x=x*(1+label.offset/norm),y=y*(1+label.offset/norm),label=label,angle=angle,hjust=hjust),
          data=co,
          color=cust("predictors.labels.color",clr),
          size=cust("predictors.labels.size",4*label.size))
      }
    }    
  }
    
  # add co-variate arrows
  if("covariates" %in% style) {
    co <- as.data.frame(cor(data$xNumeric, data$compr[,plane]))
    if(!is.null(covariates)) {
      co <- co[rownames(co)%in%covariates,]
    }
    names(co) <- c("x", "y")
    co$norm <- sqrt(co$x^2+co$y^2)
    co$label <- rownames(co)
    
    # adjust label position
    co$angle <- atan2(co$y,co$x)*180/pi
    co$hjust <- ifelse(abs(co$angle)>90,1,0)
    if(label.auto) {
      co$angle <- ifelse(abs(co$angle)>90,co$angle+180,co$angle)
    } else {
      co$angle <- 0
    }
    
    # filter according to given threshold
    if("threshold" %in% style)
      co <- co[co$norm>threshold,]
    
    if(nrow(co)==0) {
      warning("No correlation higher than threshold value ",threshold,"!")
    } else {
      clr <- cust("covariates.color","black")
      if(cust("covariates.arrows",TRUE)) {
        p <- p + geom_segment(
          aes(x=0,y=0,xend=x,yend=y),
          data=co,
          color=cust("covariates.arrows.color",clr),
          arrow=arrow(length=unit(0.02,"npc")))
      } else {
        # reset label adjust
        co$hjust<-0.5
        co$angle<-0
      }
      if(cust("covariates.labels",TRUE)) {
        p <- p + geom_text(
          aes(x=x*(1+label.offset/norm),y=y*(1+label.offset/norm),label=label,angle=angle,hjust=hjust),
          data=co,
          color=cust("covariates.labels.color",clr),
          size=4*label.size)
      }
    }
  }

  # add factors
  if("factor" %in% style) {
    bary <- aggregate(data$compr[,plane],data$xFactors[factor],mean)
    names(bary) <- c(factor,"x","y")
    if(cust("factor.points",TRUE)) {
      p <- p + geom_point(
        aes_string(x="x",y="y",color=factor),
        data=bary,
        size=cust("factor.points.size",4),
        shape=cust("factor.points.shape",13)
        )
    }
    if(cust("factor.labels",TRUE)) {
      p <- p + geom_text(
        aes_string(x="x",y="y",label=factor),
        data=bary,
        color=cust("factor.labels.color","black"),
        size=cust("factor.labels.size",4*label.size)
        )
    }
  }
  
  # return plot
  p
}

#' @title Barplot of percent of overall X variance captured by component
#' @export
#' @method barplot SCGLR
#' @S3method barplot SCGLR
#' @param height object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots optional arguments.
#' @return an object of class ggplot.
#' @seealso For barplot application see examples in \code{\link{plot.SCGLR}}.
barplot.SCGLR <- function(height, ...) {
  inertia <- data.frame(inertia=height$inertia,comp=1:length(height$inertia))
  qplot(comp, inertia, data=inertia, geom="bar", stat="identity", width=0.5) +
    scale_x_discrete(labels=names(height$inertia),limits=1:length(inertia$comp))+
    labs(x="Components", y="Inertia", title="Inertia per component\n", ...)
}

#' @title Pairwise scglr plot on components
#' @export
#' @method pairs SCGLR
#' @S3method pairs SCGLR
#' @param x object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots optionally, further arguments forwarded to \code{link{plot.SCGLR}}.
#' @param nrow number of rows of the grid layout.
#' @param ncol number of columns of the grid layout.
#' @param components vector of integers selecting components to plot (default is all components).
#' @return an object of class ggplot.
#' @seealso For pairs application see examples in \code{\link{plot.SCGLR}} 
pairs.SCGLR <- function(x, ..., nrow=NULL, ncol=NULL, components=NULL) {
  prm <- list(...)

  nr <- nrow
  nc<-ncol

  prm["x"] <- list(x)
  prm[["components"]] <- NULL
  
  # pairs of components
  if(is.null(components)) {
    ncomp<-ncol(x$compr)
  } else {
    ncomp <- components
  }
  # sanity check
  if((min(ncomp)<1) || (max(ncomp)>ncol(x$compr))) {
    stop("Invalid components for plane!")
  }
  # build pairs
  cmp_pairs <- combn(ncomp, 2, simplify=F)
  
  # build plot list
  one_plot <- function(cmp_pair) {
    do.call("plot.SCGLR", c(prm, plane=list(cmp_pair))) +
      ggtitle(NULL)
  }
  plots <- lapply(cmp_pairs, one_plot)
  
  # arrange them in a grid
  if(is.null(nc) & is.null(nr)) { 
    nr <- as.integer(sqrt(length(plots)))
  }
  do.call("arrange", c(plots, nrow=nr, ncol=nc))
}

## equivalent du par pour les ggplot2
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}