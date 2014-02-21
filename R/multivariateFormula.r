#' @title Formula construction
#' @description Construction of the multivariate scglr formula.
#' @export
#' @param namesY a vector of character containing the names of the dependent variables.
#' @param namesX a vector of character containing the names of the covariates (X) involved
#'  in the components.
#' @param namesAX a vector of character containing the names of the additional covariates.
#' @return an object of class \code{Formula}.
multivariateFormula <- function(namesY,namesX,namesAX=NULL)
{
  forml <- paste(namesY,collapse="+")
  formr <- paste(namesX,collapse="+")
  if(!is.null(namesAX)){
   formr2 <- paste(namesAX,collapse="+")
   form <- paste(paste(forml,"~",formr,sep=""),formr2,sep="|")
   formula <- as.Formula(form)
  }else{
    formula <- as.Formula(paste(forml,"~",formr,sep=""))
  }
  environment(formula) <- .GlobalEnv
  return(formula)
}