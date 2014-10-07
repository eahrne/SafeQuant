# Grapihcs.R unit tests
# 
# Author: erikahrne
###############################################################################

### INIT

source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/initTestSession.R")

### INIT END


### TEST FUNCTIONS

test.getConditionColors <- function(){
	
	cat("--- test.getConditionColors: --- \n")
	stopifnot(nrow(.getConditionColors(eset)) == 3)
	stopifnot(length(.getConditionColors(eset)[pData(eset)$condition,]) == 6)
	cat("--- test.getConditionColors: PASS ALL TEST --- \n")
}

### TEST FUNCTIONS END

### TESTS

test.getConditionColors()

if(T){
	
	pdf("/Users/erikahrne/tmp/tmp.pdf")
	# plots
	plotExpDesign(eset)
	
	# id plots 
	isDecoy <- isDecoy(fData(eset)$proteinName)
	qvals <- getIdLevelQvals(fData(eset)$idScore, isDecoy)
	idScore <- fData(eset)$idScore
	pMassError <- fData(eset)$pMassError
	
	plotScoreDistrib(idScore[!isDecoy]
			,idScore[isDecoy])
	plotIdScoreVsFDR(idScore,qvals,lwd=2)
	plotROC(qvals,breaks=seq(0,0.2,length=50)
			,main=paste("Nb. Valid identificaitons: ",sum(qvals < 0.01),"\n( FDR ",0.01,")"))
	
	# precursor mass error
	plotPrecMassErrorDistrib(pMassError,isDecoy,pMassTolWindow=c(-10,10))
	plotPrecMassErrorVsScore(pMassError, idScore, isDecoy, pMassTolWindow=c(-0.5,0.5) )
	
	# quant QC plots
	pairsAnnot(exprs(eset))
	pairsAnnot(getSignalPerCondition(eset))
	plotMSSignalDistributions(exprs(eset),col=COLORS, cex=1, cex.axis=1.5, cex.lab=1.5, ylab="binned count", xlab="AUC")
	barplotMSSignal(eset,cex.lab=1.5)
	
	##quant differential abundance related plots
	
	### plot volcanos for all case control comparisons
	plotVolcano(sqa
			,main= "created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=F
	)
	
	plotVolcano(sqa
			,main= "created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=T
	)
	
	plotVolcano(data.frame( ratio=as.vector(unlist(sqa$ratio["B"]))
					,qValue=as.vector(unlist(sqa$qValue["B"]))
					,cv=apply(sqa$cv[c("B","C")],1,max,na.rm=T))
			,caseCondition="cond B"
			,controlCondition="cond C"
			,main= "created from data.frame"
			,cex.axis=1.2
			,cex.lab=1.2
	)
	
	hClustHeatMap(eset)
	
	par(mfrow=c(2,2))
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=T,main="UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=T,main="DOWN")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=F,main="UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=F,main="DOWN")
	par(mfrow=c(1,1))
	
	plotXYDensity(exprs(eset)[,1],exprs(eset)[,2])
	plotCalibrationCurve(fit, predictorName="log10(iBAQ)")

	esetTmp <- parseProgenesisPeptideCsv(file=progenesisPeptideCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisPeptideCsvFile1))
	rtNormFactors <- getRTNormFactors(esetTmp, minFeaturesPerBin=100)
	plotRTNormSummary(rtNormFactors)
	plotRTNorm(rtNormFactors,esetTmp, samples=2)
	dev.off()
	
	
	
}
### TESTS END





















