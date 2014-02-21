### R code from vignette source 'scglrVignettes.rnw'

###################################################
### code chunk number 1: scglrVignettes.rnw:53-54
###################################################
options(width=60,prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)


###################################################
### code chunk number 2: scglrVignettes.rnw:136-137
###################################################
library("SCGLR")


###################################################
### code chunk number 3: scglrVignettes.rnw:141-142 (eval = FALSE)
###################################################
## results.scglr <- scglr(formula,data,family,K,size,offset,subset,na.action,crit)


###################################################
### code chunk number 4: scglrVignettes.rnw:150-153
###################################################
ny <- paste("y",1:2,sep="")
nx <-paste("x",1:5,sep="")
nz <- paste("z",1:3,sep="")


###################################################
### code chunk number 5: scglrVignettes.rnw:155-157
###################################################
myformula <- multivariateFormula(ny,nx,nz)
myformula


###################################################
### code chunk number 6: scglrVignettes.rnw:215-217 (eval = FALSE)
###################################################
## scglrCrossVal(formula,data,family,K,nfolds,types,size,
##                 offset,subset,na.action,crit)


###################################################
### code chunk number 7: scglrVignettes.rnw:234-237
###################################################
library("SCGLR")
data("genus")
dim(genus)


###################################################
### code chunk number 8: scglrVignettes.rnw:240-241
###################################################
names(genus)


###################################################
### code chunk number 9: scglrVignettes.rnw:245-252
###################################################
ny <- names(genus)[1:27]
sx <- which(names(genus) %in% c("geology","surface"))
nx <- names(genus)[-c(1:27,sx)] 
family <- rep("poisson",length(ny))
formula <- multivariateFormula(ny,c(nx,"I(lon*lat)"),"geology")
formula
offset <- genus$surface


###################################################
### code chunk number 10: scglrVignettes.rnw:254-257
###################################################
K <- 10
genus.cv <- scglrCrossVal(formula=formula,data=genus,family=family,
                          K=K,nfolds=5,type="mspe",offset=offset)


###################################################
### code chunk number 11: tempo
###################################################
criterion <- t(apply(genus.cv,1,function(x) x/median(x)))
criterion <- apply(criterion,2,mean) 
K.cv <- which.min(criterion)-1


###################################################
### code chunk number 12: plotCv
###################################################
plot(0:K,criterion,t="l",ylab="mspe",xlab="number of components")
Axis(side=1,at=0:K)


###################################################
### code chunk number 13: scglrVignettes.rnw:280-282
###################################################
genus.scglr<-scglr(formula=formula,data=genus,family=family,
                   K=K.cv,size=NULL,offset=offset)


###################################################
### code chunk number 14: scglrVignettes.rnw:288-289
###################################################
print(genus.scglr)


###################################################
### code chunk number 15: barplotScglr
###################################################
barplot(genus.scglr)


###################################################
### code chunk number 16: plotSimple
###################################################
plot(genus.scglr)


###################################################
### code chunk number 17: plotStyle
###################################################
plot(genus.scglr, style=c("simple","predictor","thr"), thr=0.8)


###################################################
### code chunk number 18: pairsScglr
###################################################
pairs(genus.scglr,components=c(1,3,5),ncol=2,label.size=0.5) 


