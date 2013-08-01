\name{multivariateFormula}
\alias{multivariateFormula}
\title{Formula construction}
\usage{
  multivariateFormula(namesY, namesX, namesSX = NULL)
}
\arguments{
  \item{namesY}{a vector of character containing the names
  of the dependent variables}

  \item{namesX}{a vector of character containing the names
  of the covariates (X) involved in the components}

  \item{namesSX}{a vector of character containing the names
  of the additional covariates}
}
\value{
  an object of class formula
}
\description{
  Construction of the multivariate scglr formula
}

