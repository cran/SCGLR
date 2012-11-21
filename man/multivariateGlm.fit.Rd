\name{multivariateGlm.fit}
\alias{multivariateGlm.fit}
\title{Multivariate generalized linear regression}
\usage{
  multivariateGlm.fit(Y, comp, family, offset, size)
}
\arguments{
  \item{Y}{matrix of dependent variables}

  \item{comp}{reduced scglr components obtained from the
  kComponents function}

  \item{family}{a vector of character giving the family
  distribution of each response}

  \item{offset}{used for the poisson dependent variables. A
  vector or a matrix of size: number of observations *
  number of Poisson dependent variables is expected}

  \item{size}{a matrix giving the number of trials for each
  Binomial dependent variable ncol(size) must be equal to
  the number of Binomial variables}
}
\value{
  the list, each item of which is the glm object associated
  with each response
}
\description{
  multivariateGlm is used to fit multivariate generalized
  linear models specified by a symbolic formula together
  with the distributions of the responses. This function
  performs a simple GLM fit for each dependent variable
  with the associated distribution.
}

