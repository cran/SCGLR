#' @export
#' @title Summary
#' @method summary SCGLR
#' @S3method summary SCGLR
#' @param object a SCGLR object from scglr() function
#' @param digits minimal number of significant digits, see print.default.
#' @param ... unused
summary.SCGLR <- function(object,...,digits=3) {
  print(object)
  
  dots <- list(...)
  cutoff <- ifelse(is.null(dots[["cutoff"]]),0.05,dots[["cutoff"]])
  
  rho <- as.data.frame(cor(object$xNumeric,object$compr))
  rho.pred <- as.data.frame(cor(object$lin.pred,object$compr))
  ncomp <- ncol(object$compr)
  cmp_pairs <- combn(ncomp, 2, simplify=F)

  best_plane <- function(var) {
    magni <- lapply(cmp_pairs, function(pair) sum(var[pair]^2))
    ind <- which.max(magni)
    list(bp=paste(as.character(cmp_pairs[[ind]]),collapse="/"),val=magni[ind])
  }
  tmp <- unlist(apply(rho,1,best_plane))
  tmp.pred <- unlist(apply(rho.pred,1,best_plane))
  rho <- data.frame(rho^2,best_plane=tmp[seq(1,length(tmp),2)],best_val=as.numeric(tmp[seq(2,length(tmp),2)])) 
  rho <- rho[order(rho$best_val,decreasing=TRUE),]
  cat("\nSquared correlations with numerical covariates:\n")
  print(rho,print.gap=2,digits=digits)
  rho.pred <- data.frame(rho.pred^2,best_plane=tmp.pred[seq(1,length(tmp.pred),2)],best_val=as.numeric(tmp.pred[seq(2,length(tmp.pred),2)]))
  #rho.pred <- rho.pred[order(rho.pred$best_val,decreasing=TRUE),]
  
  cat("\nSquared correlations with linear predictors:\n")  
  print(rho.pred,print.gap=2,digits=digits)
  
  coef <- sapply(object$gamma, function(x) x[,1])
  signif <- sapply(object$gamma, function(x) x[,4])
  coef[signif>=cutoff] <- NA
  colnames(coef)<-rownames(rho.pred)
  
  cat("\nGLR gamma coefficients for dependant variables:\n")
  print(coef,na.print="",digits=digits)

  invisible(list(rho,rho.pred,coef))
}