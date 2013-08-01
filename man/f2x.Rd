\name{f2x}
\alias{f2x}
\title{Regression coefficients conversion}
\usage{
  f2x(Xcr, centerx, invsqrtm, gamma, u, comp)
}
\arguments{
  \item{Xcr}{matrix of the standardized regressors}

  \item{centerx}{vector of the regressors' means}

  \item{invsqrtm}{metric matrix M^(-1/2)}

  \item{gamma}{(generalized) regression coefficients
  associated with the components}

  \item{u}{matrix of the component coordinates on Xcr
  (component = Xcr * u)}

  \item{comp}{matrix of the components}
}
\value{
  beta matrix containing the regression coefficients
  associated with the regressors X
}
\description{
  converts the regression coefficients associated with the
  components into coefficients associated with the initial
  regressors
}

