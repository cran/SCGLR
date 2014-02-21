\name{plot.SCGLR}
\alias{plot.SCGLR}
\title{SCGLR generic plot}
\usage{
\method{plot}{SCGLR}(x, ..., style = "simple", threshold = 0.8,
  plane = c(1, 2), factor = NULL, predictors = NULL, covariates = NULL,
  label.offset = 0.01, label.auto = T, label.size = 1, expand = 1)
}
\arguments{
  \item{x}{an object from SCGLR class.}

  \item{\dots}{optional arguments.}

  \item{style}{describes which plot will be drawn. Style
  "simple" : correlation plot, "covariates" (to add
  covariates arrows), "circle" (to add the unit circle),
  "observations" (to add observations), "threshold" (to
  select covariates whose sum of square correlations with
  the two components of the plane exceeds the given
  threshold), "factor" (to add one the factors),
  "predictors" (to add linear predictor arrows). Can be
  given as a vector or list of strings (eg:
  c("simple","circle")) or comma separated string (eg:
  "simple,circle"). Style elements can also be abbreviated
  (eg: "simp,cir").}

  \item{threshold}{a numeric value. All covariates whose
  sum of square correlations with the two components of the
  plane lower than this threshold will be ignored.}

  \item{plane}{a size-2 vector (or comma separated string)
  indicating which components are plotted (eg: c(1,2) or
  "1,2").}

  \item{factor}{factor to show (default to first one).}

  \item{predictors}{a vector of character to select the
  linear predictors displayed.}

  \item{covariates}{a vector of character to select the
  covariates displayed.}

  \item{label.offset}{offset by which labels should be
  moved from tip of arrows.}

  \item{label.auto}{whether or not the labels should be
  rotated according to vector angle.}

  \item{label.size}{relative size for labels (fine
  tuning).}

  \item{expand}{expand factor for windows size, for example
  to make room for clipped labels.}
}
\value{
an object of class \code{ggplot}.
}
\description{
SCGLR generic plot
}
\examples{
\dontrun{
library(SCGLR)

# load sample data
data(genus)

# get variable names from dataset
n <- names(genus)
ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
nx <- n[-grep("^gen",n)]   # X <- remaining names

# remove "geology" and "surface" from nx
# as surface is offset and we want to use geology as additional covariate
nx <-nx[!nx\%in\%c("geology","surface")]

# build multivariate formula
# we also add "lat*lon" as computed covariate
form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))

# define family
fam <- rep("poisson",length(ny))

genus.scglr <- scglr(formula=form,data = genus,family=fam, K=4,
 offset=genus$surface)

summary(genus.scglr)

barplot(genus.scglr)

plot(genus.scglr)

plot(genus.scglr,style="circle,cov,predictors,fact")

pairs(genus.scglr)

}
}

