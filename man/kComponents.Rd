\name{kComponents}
\alias{kComponents}
\title{Estimation algorithm for K components}
\usage{
  kComponents(X, Y, AX, K, family, size = NULL,
    offset = NULL, crit = list())
}
\arguments{
  \item{X}{matrix (n*p) containing the standardized
  covariates}

  \item{Y}{matrix (n*q) containing dependent variables}

  \item{AX}{matrix of additional covariates used in the
  generalized regression but not entering the linear
  combinations giving components}

  \item{K}{integer specifying the number of components}

  \item{family}{a vector of the same length as the number
  of responses containing characters identifying the
  distribution families of the dependent variables.
  "bernoulli", "binomial", "poisson" or "gaussian" are
  allowed.}

  \item{size}{matrix of size statistical units * number of
  binomial responses, giving the number of trials for
  binomial dependent variables.}

  \item{offset}{used for the poisson dependent variables. A
  vector or a matrix of size: number of observations *
  number of Poisson dependent variables is expected}

  \item{crit}{a list of maxit and tol, default is 50 and
  10e-6. If responses are bernoulli variables only, tol
  should generally be increased}
}
\value{
  a list

  \item{u}{matrix of size (number of regressors * number of
  components), contains the component-loadings, i.e. the
  coefficients of the regressors in the linear combination
  giving each component}

  \item{comp}{matrix of size (number of statistical units *
  number of components) having the components as column
  vectors}

  \item{compr}{matrix of size (number of statistical units
  * number of components) having the standardized
  components as column vectors}

  \item{ds}{the final value of the regularization degree}
}
\description{
  calculates the K components by iteratively calling
  function oneComponent
}

