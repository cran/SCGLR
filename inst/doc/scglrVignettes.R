## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
  concordance=FALSE, fig.path='scglr-',tidy=FALSE,size="small"
)

## ----echo=FALSE-------------------------------------------
options(width=60,prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)

## ----echo=FALSE-------------------------------------------
suppressPackageStartupMessages(library("SCGLR"))

## ----eval=FALSE-------------------------------------------
#  results.scglr <- scglr(formula,data,family,K,size,offset,subset,
#                         na.action,crit,method)

## ----results='hide',echo=FALSE----------------------------
ny <- paste("y",1:2,sep="")
nx <-paste("x",1:5,sep="")
nz <- paste("z",1:3,sep="")

## ----eval=TRUE--------------------------------------------
myformula <- multivariateFormula(ny,nx,nz)
myformula

## ----eval=FALSE-------------------------------------------
#  scglrCrossVal(formula,data,family,K,nfolds,types,size,
#                  offset,subset,na.action,crit,method)

## ---------------------------------------------------------
data("genus")
dim(genus)

## ---------------------------------------------------------
names(genus)

## ----eval=TRUE--------------------------------------------
ny <- names(genus)[1:27]
sx <- which(names(genus) %in% c("geology","surface"))
nx <- names(genus)[-c(1:27,sx)] 
family <- rep("poisson",length(ny))
formula <- multivariateFormula(ny,c(nx,"I(lon*lat)"),"geology")
formula
offset <- genus$surface

## ----eval=TRUE,cache=TRUE---------------------------------
K <- 15
genus.cv <- scglrCrossVal(formula=formula,data=genus,family=family,
                          K=K,nfolds=5,type="mspe",offset=offset,
                          method=methodSR(l=1, s=1/2))

## ----tempo,eval=TRUE--------------------------------------
criterion <- t(apply(genus.cv,1,function(x) x/mean(x)))
criterion.mean <- apply(criterion,2,mean) 
K.cv <- which.min(criterion.mean)-1

## ----plotCv,echo=FALSE,include=FALSE----------------------
plot(0:K,criterion.mean, type="l",
     xlab="K, number of components", ylab="Criterion (MSPE)")
Axis(side=1,at=0:K)
abline(v=K.cv,col=2)

## ----eval=TRUE,cache=TRUE---------------------------------
genus.scglr<-scglr(formula=formula,data=genus,family=family,
                   K=K.cv,size=NULL,offset=offset,
                   method=methodSR(l=1,s=1/2))

## ----eval=TRUE--------------------------------------------
print(genus.scglr)

## ----barplotScglr,include=FALSE---------------------------
barplot(genus.scglr)

## ----plotSimple,include=FALSE-----------------------------
plot(genus.scglr)

## ----plotStyle,include=FALSE------------------------------
plot(genus.scglr, threshold=0.8, predictors=TRUE)

## ----pairsScglr,include=FALSE,eval=TRUE-------------------
pairs(genus.scglr,components=c(1,3,5),ncol=2,label.size=0.5) 

