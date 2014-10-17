# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

myRatioBoxplot <- function(eset, sampleNb=nrow(eset),xlab="Reference Ratio",ylab="TMT Ratio", ...){
	
	ratios <-  2^getRatios(eset)
	sampleNb <- min(c(nrow(eset),sampleNb))
	sampleIndices <- sample(1:nrow(ratios),sampleNb,replace=F)
	par(mar=c(5.1,5.1,4.1,2.1))
	boxplot(ratios[sampleIndices,c(2,4,1,3)],log="y", ylim=c(0.1,25),names=c("0.2:1","0.5:1","1.5:1","10:1"),lwd=2, col=COLORS[2:5],..., cex.lab=1.5,cex.axis=1.5, xlab=xlab,ylab=ylab)
	abline(h=c(0.2,0.5,1.5,10), col=COLORS[2:5],lty=2)
	abline(h=1,lty=1)
	par(mfrow=c(1,1))
}

source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")


load("/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda")


sqaSpectrum <- safeQuantAnalysis(esetTMT6Spectrum, method="global")


myRatioBoxplot(sqaSpectrum$eset[!fData(sqaSpectrum$eset)$isNormAnchor,])

par(mfrow=c(1,1))
boxplot((2^sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,])-1)






pdf("/Users/erikahrne/tmp/tmp.pdf")

par(mfrow=c(1,2))

### ratios 1.5,0.2,10,0.5
medianTMTRatio <-  apply(sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,] ,2,median, na.rm=T )
refRatio <-  log2(c(1.5,0.2,10,0.5))
plot(refRatio,medianTMTRatio)
abline(lm(medianTMTRatio ~ refRatio))
abline(coef=c(0,1), lty=2)

expDesignPosRatios <- data.frame(condition=paste("cond",c(1,2,3,1,4,5),sep="_"),isControl=rep(F,6) )
expDesignPosRatios$isControl[c(3)] <- T


esetPosRatios <- createExpressionDataset(expressionMatrix=exprs(esetTMT6Spectrum),expDesign=expDesignPosRatios,featureAnnotations=fData(esetTMT6Spectrum))

posRatios <- apply(getRatios(esetPosRatios,log2=T)[!fData(sqaSpectrum$eset)$isNormAnchor,] ,2,median, na.rm=T )


refRatio <-  log2(c(5,7.5,50,2.5))
plot(refRatio,posRatios)
abline(lm(posRatios ~ refRatio))
abline(coef=c(0,1), lty=2)

dev.off()

print("DONE")





