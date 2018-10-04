# Grapihcs.R unit tests
#
# Author: ahrnee-adm
###############################################################################

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
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

testPlotNbIdentificationsVsRT <- function(){

	file <- progenesisPeptideMeasurementFile1
	#file <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements5_MRuegg.csv"

	cat("--- testPlotNbIdentificationsVsRT: --- \n")
	esetTmp <- parseProgenesisPeptideMeasurementCsv(file,expDesign= getExpDesignProgenesisCsv(file))
	plotNbIdentificationsVsRT(esetTmp)
	cat("--- testPlotNbIdentificationsVsRT: PASS ALL TEST --- \n")

}

if(T){

	tmpFile <- paste(tempdir(),"/tmp.pdf",collapse="",sep="")

	pdf(tmpFile)
	# plots
	plotExpDesign(eset)

	# id plots
	isDec <- isDecoy(fData(eset)$proteinName)
	qvals <- getIdLevelQvals(fData(eset)$idScore, isDec)
	idScore <- fData(eset)$idScore
	pMassError <- fData(eset)$pMassError

	plotScoreDistrib(idScore[!isDec]
			,idScore[isDec], main="plotScoreDistrib")
	plotIdScoreVsFDR(idScore,qvals,lwd=2, main="plotIdScoreVsFDR")
	plotROC(qvals,breaks=seq(0,0.2,length=50)
			,main=paste("plotROC Nb. Valid identificaitons: ",sum(qvals < 0.01),"\n( FDR ",0.01,")"))

	# precursor mass error
	plotPrecMassErrorDistrib(eset,pMassTolWindow=c(-10,10),main="plotPrecMassErrorDistrib")
	plotPrecMassErrorVsScore(eset, pMassTolWindow=c(-0.5,0.5),main="plotPrecMassErrorVsScore")

	# quant QC plots
	pairsAnnot(exprs(eset), main="pairsAnnot")
	pairsAnnot(getSignalPerCondition(eset), main="pairsAnnot")

	plotMSSignalDistributions(exprs(eset),col=COLORS, cex=1, cex.axis=1.5, cex.lab=1.5, ylab="binned count", xlab="AUC", main="plotMSSignalDistributions")
	barplotMSSignal(eset,cex.lab=1.5,main="barplotMSSignal")

	### TMT calibration mix
#	.plotTMTRatioVsRefRatio(rollUp(esetCalibMixPair , featureDataColumnName="peptide"), cex.lab=1.5, cex.axis=1.5)
#
	##quant differential abundance related plots

	### plot volcanos for all case control comparisons
	plotVolcano(sqa
			,main= "plotVolcano created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=F
	)

	plotVolcano(sqa
			,main= "plotVolcano created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=T
	)

	plotVolcano(data.frame( ratio=as.vector(unlist(sqa$ratio["B"]))
					,qValue=as.vector(unlist(sqa$qValue["B"]))
					,cv=apply(sqa$cv[c("B","C")],1,max,na.rm=T))
			,caseCondition="cond B"
			,controlCondition="cond C"
			,main= "plotVolcano created from data.frame"
			,cex.axis=1.2
			,cex.lab=1.2
	)

	plotVolcano(safeQuantAnalysis(esetPaired)
			,main= "plotVolcano created from safeQuantAnalysis object (esetPaired)"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=F
	)

	hClustHeatMap(eset)
	hClustHeatMap(eset, dendogram="both")

	par(mfrow=c(2,2))
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=T,main="plotNbValidDeFeaturesPerFDR UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=T,main="plotNbValidDeFeaturesPerFDR DOWN")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=F,main="plotNbValidDeFeaturesPerFDR UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=F,main="plotNbValidDeFeaturesPerFDR DOWN")
	par(mfrow=c(1,1))

	plotXYDensity(exprs(eset)[,1],exprs(eset)[,2], main="plotXYDensity")
	plotAbsEstCalibrationCurve(absEstSimDataFit, predictorName="log10(iBAQ)", main="plotCalibrationCurve")

	esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	rtNormFactors <- getRTNormFactors(esetTmp, minFeaturesPerBin=100)
	plotRTNormSummary(esetTmp, main="plotRTNormSummary")
	plotRTNorm(rtNormFactors,esetTmp, samples=2, main="plotRTNorm")

	missinValueBarplot(eset)
	plotQValueVsPValue(sqa, lim=c(0,0.5))
	.correlationPlot(exprs(eset))

	maPlotSQ(eset)

	eset <-  parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementCsvFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementCsvFile1))
	cysteinFreqBarplot(fData(eset)$peptide)

	graphics.off()

	cat("CREATED FILE", tmpFile, "\n")

}
### TESTS END

# library(SafeQuant)
# library(stringr)
#
# ecoliDb =  read.fasta("/Volumes/pcf01$/Schmidt_Group/Databases/Latest_UniProt/uniprot-ecoli.fasta",seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
# humanDb =  read.fasta("/Volumes/pcf01$/Schmidt_Group/Databases/Latest_UniProt/uniprot-Human_301014.fasta",seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
# mouseDb =  read.fasta("/Volumes/pcf01$/Schmidt_Group/Databases/Latest_UniProt/uniprot-Mouse_301014.fasta",seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
# yeastDb =  read.fasta("/Volumes/pcf01$/Schmidt_Group/Databases/Latest_UniProt/uniprot-yeast_301014.fasta",seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
#
#
# ecoliPeptides = lapply(ecoliDb, getPeptides ) %>% unlist %>% unique
# ecoliPeptidesLength = nchar(ecoliPeptides)
# ecoliPeptides = subset(ecoliPeptides ,ecoliPeptidesLength > 6 & ecoliPeptidesLength < 20  )
# ecoliCFreq = str_detect(ecoliPeptides,"C") %>% sum / length(ecoliPeptides)
#
# humanPeptides = lapply(humanDb, getPeptides ) %>% unlist %>% unique
# humanPeptidesLength = nchar(humanPeptides)
# humanPeptides = subset(humanPeptides ,humanPeptidesLength > 6 & humanPeptidesLength < 20  )
# humanCFreq = str_detect(humanPeptides,"C") %>% sum / length(humanPeptides)
#
# mousePeptides = lapply(mouseDb, getPeptides ) %>% unlist %>% unique
# mousePeptidesLength = nchar(mousePeptides)
# mousePeptides = subset(mousePeptides ,mousePeptidesLength > 6 & mousePeptidesLength < 20  )
# mouseCFreq = str_detect(mousePeptides,"C") %>% sum / length(mousePeptides)
#
# yeastPeptides = lapply(yeastDb, getPeptides ) %>% unlist %>% unique
# yeastPeptidesLength = nchar(yeastPeptides)
# yeastPeptides = subset(yeastPeptides ,yeastPeptidesLength > 6 & yeastPeptidesLength < 20  )
# yeastCFreq = str_detect(yeastPeptides,"C") %>% sum / length(yeastPeptides)

#getNbMisCleavages(fData(esetNorm)$peptide)




# file = "/Users/ahrnee-adm/dev/R/workspace/SafeQuantTestData/Progenesis/v2/peptide_measurements2.csv"
#
# eset = parseProgenesisPeptideMeasurementCsv(file,expDesign=getExpDesignProgenesisCsv(file) )
#
# eset = sqNormalize(eset)
# eset = sqImpute(eset)
# esetProt = rollUp(eset)
#
# #eset <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisProteinCsvFile1))
#
#
# hClustHeatMap(eset)



