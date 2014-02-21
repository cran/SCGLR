\name{critConvergence}
\alias{critConvergence}
\title{Auxiliary for controlling SCGLR fitting}
\usage{
critConvergence(tol = 1e-06, maxit = 50)
}
\arguments{
  \item{tol}{positive convergence threshold.}

  \item{maxit}{integer, maximum number of iterations.}
}
\value{
a list containing elements named as the arguments.
}
\description{
Auxiliary function for \code{scglr} fitting used to
construct a convergence control argument.
}

