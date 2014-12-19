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

#expDesignPosRatios <- data.frame(condition=paste("cond",c(1,2,3,1,4,5),sep="_"),isControl=rep(F,6) )
#expDesignPosRatios$isControl[c(3)] <- T
#
#
#esetPosRatios <- createExpressionDataset(expressionMatrix=exprs(esetTMT6Spectrum),expDesign=expDesignPosRatios,featureAnnotations=fData(esetTMT6Spectrum))
#posRatios <- apply(getRatios(esetPosRatios,log2=T)[!fData(sqaSpectrum$eset)$isNormAnchor,] ,2,median, na.rm=T )
#
#refRatio <-  log2(c(5,7.5,50,2.5))
#plot(refRatio,posRatios)
#abline(lm(posRatios ~ refRatio))
#abline(coef=c(0,1), lty=2)


### INTERFERENCE

tmtMedianInt <- getSignalPerCondition(sqaSpectrum$eset)
ISBARHE <- !fData(sqaSpectrum$eset)$isNormAnchor

# calcualate interference level
interferenceLevels <-   data.frame((tmtMedianInt[,2]-(tmtMedianInt[,1]*1.5) )/(1- 1.5)
		,(tmtMedianInt[,3]-(tmtMedianInt[,1]*0.2) )/(1- 0.2)
		,(tmtMedianInt[,4]-(tmtMedianInt[,1]*10) )/(1- 10)
		,(tmtMedianInt[,5]-(tmtMedianInt[,1]*0.5) )/(1- 0.5)

)
names(interferenceLevels) <- c("r1.5","r0.2","r10","r0.5")


boxplot(log(interferenceLevels))

d <-  tmtMedianInt[,3] - interferenceLevels
d[ISBARHE,2]


### get median interference level
### weigh in accordance with BARHE:HUMAN ratio and normalization factors
interferenceLevelMedian <- apply(interferenceLevels,1,median)
#interferenceLevelMedian <- interferenceLevels[,2]

### adj. interference levels to never be higher than minimal TMT int
nonOk <- (interferenceLevelMedian > apply(exprs(sqaSpectrum$eset),1,min) ) | is.na(interferenceLevelMedian) | (interferenceLevelMedian < 0)
interferenceLevelMedian[ nonOk ]  <- apply(exprs(sqaSpectrum$eset),1,min) [ nonOk ]
#interferenceLevelMedian  <- interferenceLevelMedian[ !nonOk ] 

### get relative interference levels 
relativeInterference <- (interferenceLevelMedian / exprs(sqaSpectrum$eset))




plotXYDensity( log(fData(sqaSpectrum$eset)$ms1Area[ISBARHE]) ,relativeInterference[ISBARHE,1] )

plotXYDensity( log(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]) ,relativeInterference[ISBARHE,1] )

plotXYDensity( log(apply(exprs(sqaSpectrum$eset),1,mean)[ISBARHE]) ,relativeInterference[ISBARHE,1] )


plotXYDensity( as.vector(unlist(log(exprs(sqaSpectrum$eset[ISBARHE,])))) , as.vector(unlist(relativeInterference[ISBARHE,])) )
plotXYDensity( as.vector(unlist(log(exprs(sqaSpectrum$eset[ISBARHE,]) * fData(sqaSpectrum$eset)$injectionTime[ISBARHE] ))) , as.vector(unlist(relativeInterference[ISBARHE,])) )



plotXYDensity( log2(fData(sqaSpectrum$eset)$ms1Area[ISBARHE] ) ,log2(apply(exprs(sqaSpectrum$eset),1,sum)[ISBARHE]) )


plotXYDensity( log2(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]) ,log2(apply(exprs(sqaSpectrum$eset),1,sum)[ISBARHE]) )


ok <- (interferenceLevelMedian != apply(exprs(sqaSpectrum$eset),1,min))

plotXYDensity( log2( apply(exprs(sqaSpectrum$eset)[ISBARHE & ok,],1,sum)/ fData(sqaSpectrum$eset)$injectionTime[ISBARHE & ok] ) , log2(interferenceLevelMedian[ISBARHE & ok]))
plotXYDensity( log2(apply(exprs(sqaSpectrum$eset)[ISBARHE & ok,],1,sum) / fData(sqaSpectrum$eset)$injectionTime[ISBARHE & ok] ) , log2(interferenceLevelMedian[ISBARHE & ok]))


plotXYDensity(log2(fData(sqaSpectrum$eset)$ms1Int[ISBARHE & ok] / fData(sqaSpectrum$eset)$injectionTime[ISBARHE & ok]   ) , log2(interferenceLevelMedian[ISBARHE & ok]))


y <- 1/sqrt(rep(interferenceLevelMedian[ISBARHE & ok],6))
x <- as.vector(unlist(log2(exprs(sqaSpectrum$eset[ISBARHE & ok,]))))

y <- log(rep(interferenceLevelMedian[ISBARHE ],6))
x <- as.vector(unlist(exprs(sqaSpectrum$eset[ISBARHE,])))

plotXYDensity(x  ,y  )

library(MASS)
fit <- rlm(y ~x, weights=1/sqrt(x) ) 


abline(fit)

pred <- predict(fit)

plotXYDensity(y[as.numeric(names(pred))],pred)

qqnorm(fit$resid,xlim=c(-3,3),ylim=c(-3,3))
abline(coef=c(0,1))

boxplot(log2(exprs(sqaSpectrum$eset[ISBARHE & ok,])))
boxplot(relativeInterference[ISBARHE & ok ,])
boxplot(log2(interferenceLevels[ISBARHE & ok ,]))

plotXYDensity( log2( as.vector(unlist(exprs(sqaSpectrum$eset[ISBARHE & ok,]))) / rep(fData(sqaSpectrum$eset)$injectionTime[ISBARHE & ok],6)  ) , log2(rep(interferenceLevelMedian[ISBARHE & ok],6)) )


boxplot( relativeInterference[ISBARHE,1] ~ round(log2(fData(sqaSpectrum$eset)$ms1Area[ISBARHE])) )
boxplot( relativeInterference[ISBARHE,1] ~ round(log2(fData(sqaSpectrum$eset)$ms1Int[ISBARHE])) )
boxplot( relativeInterference[ISBARHE,6] ~ round(log2(apply(exprs(sqaSpectrum$eset),1,max)[ISBARHE])) )


plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,sum)[ISBARHE]), relativeInterference[ISBARHE,6] )

plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,max)[ISBARHE]), relativeInterference[ISBARHE,6] )

plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,min)[ISBARHE]), relativeInterference[ISBARHE,6] )


plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,min)[ISBARHE]), log2(interferenceLevelMedian[ISBARHE] ))
plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,max)[ISBARHE]), log2(interferenceLevelMedian[ISBARHE] ))
plotXYDensity( log2(apply(exprs(sqaSpectrum$eset),1,sum)[ISBARHE]), log2(interferenceLevelMedian[ISBARHE] ))


plotXYDensity(log2(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]/fData(sqaSpectrum$eset)$injectionTime[ISBARHE]) , relativeInterference[ISBARHE,1]  )


boxplot(log2(exprs(sqaSpectrum$eset))[ISBARHE,])


boxplot( relativeInterference[ISBARHE,6] ~ round(log2(fData(sqaSpectrum$eset)$ms1Int[ISBARHE])) )
boxplot( relativeInterference[ISBARHE,6] ~ round(log2(apply(exprs(sqaSpectrum$eset),1,min)[ISBARHE])) )
boxplot( relativeInterference[ISBARHE,6] ~ round(log2(apply(exprs(sqaSpectrum$eset),1,max)[ISBARHE])) )




boxplot( relativeInterference[ISBARHE,3] ~ I((10*(1+round(fData(sqaSpectrum$eset)$injectionTime[ISBARHE]/10)))))

boxplot( relativeInterference[ISBARHE,5] ~ I((10*(1+round(fData(sqaSpectrum$eset)$interference[ISBARHE]/10)))))


plotXYDensity( log(fData(sqaSpectrum$eset)$ms1Area[ISBARHE]) ,log(interferenceLevelMedian[ISBARHE]) )

plotXYDensity( log(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]) ,log(interferenceLevelMedian[ISBARHE]) )

plotXYDensity(log(interferenceLevelMedian[ISBARHE]),  log(fData(sqaSpectrum$eset)$ms1Area[ISBARHE]) )

plotXYDensity(log(interferenceLevelMedian[ISBARHE]),  log(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]) )

plotXYDensity(log(fData(sqaSpectrum$eset)$ms1Int[ISBARHE]), relativeInterference[ISBARHE,5] )

library(MASS)

ok <- !(interferenceLevelMedian == apply(exprs(sqaSpectrum$eset),1,min))
i <- interferenceLevelMedian[ISBARHE & ok]
a <- apply(exprs(sqaSpectrum$eset),1,min)[ISBARHE & ok]
#b <- fData(sqaSpectrum$eset)$ms1Int[ISBARHE & ok]
b <- apply(exprs(sqaSpectrum$eset),1,max)[ISBARHE & ok]
#b <- 

#fit <- lm( i ~ a*b, weights=1/(i^2) )
fit <- rlm( i ~ a*b, weights=1/(i^2))


summary(fit)

pred <- predict(fit)

plotXYDensity(log10(i[as.numeric(names(pred))]),log10(pred))

plotXYDensity(log10(i[as.numeric(names(pred))]),log10(pred))


plotXYDensity(log10(a),log10(i))
plotXYDensity(log10(b),log10(i))


boxplot((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,1])) ~ fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$nbRolledFeatures)



plotXYDensity(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,3]))) , fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$nbRolledFeatures, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) , fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$nbRolledFeatures, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) , fData(sqaSpectrum$eset[ISBARHE,])$ms1Int, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) , fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$ms1Int, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) , fData(sqaSpectrum$eset[ISBARHE,])$interference, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) , fData(sqaSpectrum$eset[ISBARHE,])$injectionTime, use="complete.obs")


ok <- log(fData(sqaSpectrum$eset)$ms1Int) > 0

cor(unlist((interferenceLevelMedian[ISBARHE & ok] / exprs(sqaSpectrum$eset[ISBARHE & ok ,]))) , fData(sqaSpectrum$eset[ISBARHE & ok,])$ms1Int/ fData(sqaSpectrum$eset[ISBARHE& ok,])$injectionTime, use="complete.obs")

plotXYDensity(log((interferenceLevelMedian[ISBARHE & ok] )) , log(fData(sqaSpectrum$eset[ISBARHE & ok,])$ms1Int / fData(sqaSpectrum$eset[ISBARHE & ok,])$injectionTime), use="complete.obs")

plotXYDensity(log((interferenceLevelMedian[ISBARHE & ok] )) ,log(fData(sqaSpectrum$eset[ISBARHE & ok,])$ms1Int), use="complete.obs")

plotXYDensity(log((interferenceLevelMedian[ISBARHE & ok] )) , fData(sqaSpectrum$eset[ISBARHE & ok,])$injectionTime, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE & ok] / exprs(sqaSpectrum$eset[ISBARHE & ok ,]))) , fData(sqaSpectrum$eset[ISBARHE & ok,])$ms1Int * fData(sqaSpectrum$eset[ISBARHE& ok,])$injectionTime, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE & ok] / exprs(sqaSpectrum$eset[ISBARHE & ok ,]))) , fData(sqaSpectrum$eset[ISBARHE & ok,])$ms1Int * fData(sqaSpectrum$eset[ISBARHE& ok,])$injectionTime, use="complete.obs")

cor(unlist((interferenceLevelMedian[ISBARHE & ok] / exprs(sqaSpectrum$eset[ISBARHE & ok,]))) , (fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$ms1Int / fData(esetTMT6Protein[as.character(fData(sqaSpectrum$eset[ISBARHE,])$proteinName),])$interference)[ok], use="pairwise.complete.obs")^2

cor(unlist((interferenceLevelMedian[ISBARHE] / exprs(sqaSpectrum$eset[ISBARHE,]))) ,apply(exprs(sqaSpectrum$eset)[ISBARHE,],1,sum), use="complete.obs")^2




plotXYDensity(log2(fData(esetTMT6Protein)$ms1Int), log2(fData(esetTMT6Protein)$nbRolledFeatures))
plotXYDensity(log2(fData(esetTMT6Protein)$ms1Int)[!fData(esetTMT6Protein)$isNormAnchor], log2(fData(esetTMT6Protein)$nbRolledFeatures)[!fData(esetTMT6Protein)$isNormAnchor])


plotXYDensity(log2(apply(exprs(esetTMT6Protein),1,sum)/fData(esetTMT6Protein)$nbRolledFeatures), log2(fData(esetTMT6Protein)$nbRolledFeatures))

plotXYDensity( log2(fData(esetTMT6Protein)$nbRolledFeatures)[!fData(esetTMT6Protein)$isNormAnchor] , log2(apply(exprs(esetTMT6Protein),1,sum)[!fData(esetTMT6Protein)$isNormAnchor]/fData(esetTMT6Protein)$nbRolledFeatures)[!fData(esetTMT6Protein)$isNormAnchor])


plotXYDensity(log2(apply(exprs(esetTMT6Spectrum),1,sum))[!fData(esetTMT6Spectrum)$isNormAnchor], log2(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName)[!fData(esetTMT6Spectrum)$isNormAnchor],])$nbRolledFeatures) )

plotXYDensity(log2(fData(esetTMT6Spectrum)$ms1Int), log2(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures) )

plotXYDensity(log2(apply(exprs(esetTMT6Spectrum),1,sum)),log2(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$ms1Int))

plotXYDensity(log2(apply(exprs(esetTMT6Spectrum),1,min)/fData(esetTMT6Spectrum)$injectionTime), log2(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures) )
plotXYDensity(log2(apply(exprs(esetTMT6Spectrum),1,min)), log2(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures) )


boxplot(log2(exprs(esetTMT6Spectrum))[!fData(sqaSpectrum$eset)$isNormAnchor,])




plotXYDensity(log10(fData(esetTMT6Spectrum)$injectionTime), log10(fData(esetTMT6Spectrum)$ms1Int)  )

plotXYDensity(log10(fData(esetTMT6Spectrum)$injectionTime), log10(log2(apply(exprs(esetTMT6Spectrum),1,min)))  )

plotXYDensity(log10(fData(esetTMT6Spectrum)$injectionTime),log10(log2(apply(exprs(esetTMT6Spectrum),1,sum)))  )


plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,sum)))  )
plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,min)))  )



plotXYDensity(fData(esetTMT6Spectrum)$injectionTime, fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures )


#### TO BE KEPT
plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,sum)))  )
plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,sum)) / fData(esetTMT6Spectrum)$injectionTime )  )

plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,min)))  )
plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(log2(apply(exprs(esetTMT6Spectrum),1,min)) / fData(esetTMT6Spectrum)$injectionTime )  )
#### TO BE KEPT END


plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )
plotXYDensity(log10(fData(esetTMT6Spectrum)$ms1Int/ fData(esetTMT6Spectrum)$injectionTime),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )

plotXYDensity(log10(apply(exprs(esetTMT6Spectrum),1,min)),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )
plotXYDensity(log10(apply(exprs(esetTMT6Spectrum),1,min)/ fData(esetTMT6Spectrum)$injectionTime),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )

plotXYDensity(log10(apply(exprs(esetTMT6Spectrum),1,sum)),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )
plotXYDensity(log10(apply(exprs(esetTMT6Spectrum),1,sum)/ fData(esetTMT6Spectrum)$injectionTime),log10(fData(esetTMT6Protein[as.character(fData(esetTMT6Spectrum)$proteinName),])$nbRolledFeatures)  )


plotXYDensity(log10(fData(esetTMT6Spectrum)$injectionTime), log10(fData(esetTMT6Spectrum)$interference))
names(fData(esetTMT6Spectrum))
dev.off()

print("DONE")





