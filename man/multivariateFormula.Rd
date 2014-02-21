\name{multivariateFormula}
\alias{multivariateFormula}
\title{Formula construction}
\usage{
multivariateFormula(namesY, namesX, namesAX = NULL)
}
\arguments{
  \item{namesY}{a vector of character containing the names
  of the dependent variables.}

  \item{namesX}{a vector of character containing the names
  of the covariates (X) involved in the components.}

  \item{namesAX}{a vector of character containing the names
  of the additional covariates.}
}
\value{
an object of class \code{Formula}.
}
\description{
Construction of the multivariate scglr formula.
}

