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


medianTMTFC <-  apply((2^abs(sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,])) ,2,median, na.rm=T )-1

FC <-  c(1.5,5,10,2)-1

k <- c(0.4,0.1,0.4,0.3)

plot(medianTMTFC,FC)

pData(sqaSpectrum$eset)

sqaSpectrumMedianInt <- getSignalPerCondition(sqaSpectrum$eset)

par(mfrow=c(2,2))
summary(plotXYDensity(log2(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"]))

par(mfrow=c(2,2))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_1"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]))
summary(plotXYDensity(log(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_1"]),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"]))


par(mfrow=c(2,2))
summary(plotXYDensity(log2(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum))  ,sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]))

summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"]))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"]))
par(mfrow=c(1,1))



refRatios <- log2(data.frame( rep(1.5,nrow(esetTMT6Spectrum))
								,rep(0.2,nrow(esetTMT6Spectrum))
								,rep(10,nrow(esetTMT6Spectrum))
								,rep(0.5,nrow(esetTMT6Spectrum))
		
		))

relativeCompression <- sqaSpectrum$ratio

relativeCompression[,1] <-  refRatios[,1] / sqaSpectrum$ratio[,1]
relativeCompression[,2] <-  refRatios[,2] / sqaSpectrum$ratio[,2]
relativeCompression[,3] <-  refRatios[,3] / sqaSpectrum$ratio[,3]
relativeCompression[,4] <-  refRatios[,4] / sqaSpectrum$ratio[,4]

relativeCompression[relativeCompression < 0] <- NA

names(relativeCompression) <- names(sqaSpectrum$ratio)

pdf("/Users/erikahrne/tmp/tmp.pdf")

par(mfrow=c(2,2))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"], ylim=c(-10,10) ))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"], ylim=c(-10,10) ))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"], ylim=c(-10,10) ))
summary(plotXYDensity(log(apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,],1,sum)),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"], ylim=c(-10,10) ))



par(mfrow=c(2,2))
summary(plotXYDensity(log(fData(sqaSpectrum$eset)$ms1Int[!fData(sqaSpectrum$eset)$isNormAnchor]),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_2"]))
summary(plotXYDensity(log(fData(sqaSpectrum$eset)$ms1Int[!fData(sqaSpectrum$eset)$isNormAnchor]),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_3"]))
summary(plotXYDensity(log(fData(sqaSpectrum$eset)$ms1Int[!fData(sqaSpectrum$eset)$isNormAnchor]),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_4"]))
summary(plotXYDensity(log(fData(sqaSpectrum$eset)$ms1Int[!fData(sqaSpectrum$eset)$isNormAnchor]),relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,"cond_5"]))

dev.off()

par(mfrow=c(2,2))

boxplot(abs(sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,]))
boxplot(log2(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,2:5]))
plot(log2(c(1.5,0.2,10,0.5)), log2(c(1.5,0.2,10,0.5)) / apply(sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,],2,median,na.rm=T) ,ylab="Correction Factor"  )

boxplot(relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,], ylim=c(0,20))

par(mfrow=c(2,2))
plot(log2(c(1.5,0.2,10,0.5)), apply(sqaSpectrum$ratio[!fData(sqaSpectrum$eset)$isNormAnchor,],2,median,na.rm=T),ylab="TMT Ratio", xlab="Ratio")
abline(coef=c(0,1))
plot(log2(c(1.5,0.2,10,0.5)), apply(relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,],2,median,na.rm=T),ylab="Correction Factor", xlab="Ratio")
plot(log2( apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,2:5],2,median,na.rm=T) ), apply(relativeCompression[!fData(sqaSpectrum$eset)$isNormAnchor,],2,median,na.rm=T),ylab="Correction Factor", xlab="tmt median int")
plot(log2(c(1.5,0.2,10,0.5)),log2( apply(sqaSpectrumMedianInt[!fData(sqaSpectrum$eset)$isNormAnchor,2:5],2,median,na.rm=T) ), xlab="Ratio", ylab="tmt median int")

print("DONE")





