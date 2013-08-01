\name{oneComponent}
\alias{oneComponent}
\title{Function that calculates the current scglr-component}
\usage{
  oneComponent(X, Y, AX, family, size = NULL,
    offset = NULL, ds, crit)
}
\arguments{
  \item{X}{a matrix (n*p) containing covariates
  standardized through M^(-1/2)}

  \item{Y}{a matrix (n*q) containing dependent variables}

  \item{AX}{a matrix containing additional covariates}

  \item{family}{a vector of charater of length q specifying
  the distributions of the responses. Bernoulli, binomial,
  Poisson and gaussian are allowed.}

  \item{size}{specifies the numbers of trials of the
  binomial responses.  A (n*qb) matrix is expected for qb
  binomial variables.}

  \item{offset}{used for the poisson dependent variables. A
  vector or a matrix of size: number of observations *
  number of Poisson dependent variables is expected}

  \item{ds}{an integer specifying the degree to which
  structural strength in X is taken into account when
  calculating the component}

  \item{crit}{a list of maxit and tol, default is 50 and
  10e-6. If responses are bernoulli variables only, tol
  should generally be increased.}
}
\value{
  the unit vector of loadings associated with the current
  component, i.e. the coefficients of the regressors in the
  linear combination giving each component
}
\description{
  Function that calculates the current scglr-component
}

