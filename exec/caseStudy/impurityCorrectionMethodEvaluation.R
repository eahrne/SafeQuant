
# MSnBase
# http://www.bioconductor.org/packages/release/bioc/manuals/MSnbase/man/MSnbase.pdf
# Purity correction here is applied using solve from the base package using the purity correction
# values as coefficient of the linear system and the reporter quantities as the right-hand side of the
# linear system. NA values are ignored and negative intensities after correction are also set to NA.
# A more elaborated purity correction method is described in Shadforth et al., i-Tracker: for quantitative
# proteomics using iTRAQ. BMC Genomics. 2005 Oct 20;6:145. (PMID 16242023).
# 
# Author: erikahrne
###############################################################################

#@TEMP
source("/Users/erikahrne/dev/R/workspace/TMTRatioCorrection/Functions.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")
library(gplots) # volcano plot
library(data.table)

### PARAMS

#scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Raw Data Report for TMT-ratio-bart-only-200415.xls"
scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_6plex_bart_only_220515/Scaffold/Raw Data Report for TMT-ratio-bart-only-220515.xls"
rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/impurityCorrectionMethodEvaluation.rData"

### PARAMS END 


if(F){
	#expDesignTMTSixPlex <- data.frame(condition=paste("cond",c(1,2,3,1,4,5),sep="_"),isControl=rep(F,6) )
	expDesignTMTSixPlex <- data.frame(condition=paste("cond",c(1,2,3,4,5,1),sep="_"),isControl=rep(F,6) )
	
	expDesignTMTSixPlex$isControl[c(1,6)] <- T
	#expDesignTMTSixPlex$isControl[c(1,4)] <- T
	
	eset <- parseScaffoldRawFile(scaffoldRawDataFile, expDesign=expDesignTMTSixPlex, isPurityCorrect=F)
	#filter keep only BARHE
	eset <- eset[grepl("BARHE",fData(eset)$proteinName),] 
	save(file=rDataFile,eset)
}else{
	load(rDataFile)
}




### standard purity correction, set invalid to zero
esetSolveAllZero <-  eset
#esetSolveAllZero <- eset[apply(esetSolveAllZero,1,function(t){ return(sum(is.na(t)) == 0 )}),]
exprs(esetSolveAllZero) <- purityCorrectTMT(exprs(esetSolveAllZero),getImpuritiesMatrix(6), invalidReplace = "allZero")
esetProtSolveAllZero <- rollUp(esetSolveAllZero,featureDataColumnName= c("proteinName")) 

### standard purity correction, set invalid to NA
esetSolveAllNA <-  eset
exprs(esetSolveAllNA) <- purityCorrectTMT(exprs(esetSolveAllNA),getImpuritiesMatrix(6), invalidReplace = "allNA")
esetProtSolveAllNA <- rollUp(esetSolveAllNA,featureDataColumnName= c("proteinName")) 

### standard purity correction, replace invalid with original
esetSolveAllOrg <-  eset
exprs(esetSolveAllOrg) <- purityCorrectTMT(exprs(esetSolveAllOrg),getImpuritiesMatrix(6), invalidReplace = "allOrg")
esetProtSolveAllOrg <- rollUp(esetSolveAllOrg,featureDataColumnName= c("proteinName")) 

### standard purity correction, keep invalid
esetSolveKeepInvalid <-  eset
exprs(esetSolveKeepInvalid) <- purityCorrectTMT(exprs(esetSolveKeepInvalid),getImpuritiesMatrix(6), invalidReplace = "")
esetProtSolveKeepInvalid <- rollUp(esetSolveKeepInvalid,featureDataColumnName= c("proteinName")) 

### Simply subtract channel 5, keep invalid
esetSubtractCh5OnlyKeepInvalid <-  eset
exprs(esetSubtractCh5OnlyKeepInvalid)[,5]  <- exprs(esetSubtractCh5OnlyKeepInvalid)[,5] - exprs(esetSubtractCh5OnlyKeepInvalid)[,4]*0.041 - exprs(esetSubtractCh5OnlyKeepInvalid)[,6]*0.033
#exprs(esetSubtractCh5OnlyKeepInvalid)[,5]  <- exprs(eset)[,5] - exprs(esetSubtractCh5OnlyKeepInvalid)[,4]*0.041
esetSubtractCh5OnlyKeepInvalid <- esetSubtractCh5OnlyKeepInvalid[apply(esetSubtractCh5OnlyKeepInvalid,1,function(t){ return(sum(is.na(t)) == 0 )}),]
esetProtSubtractCh5OnlyKeepInvalid <- rollUp(esetSubtractCh5OnlyKeepInvalid,featureDataColumnName= c("proteinName")) 

### standard purity correction, keep invalid
esetSolveCh5OnlyKeepInvalid <-  eset
exprs(esetSolveCh5OnlyKeepInvalid) <- purityCorrectTMT(exprs(esetSolveCh5OnlyKeepInvalid),getImpuritiesMatrix(-6), invalidReplace = "")
esetProtSolveCh5OnlyKeepInvalid <- rollUp(esetSolveCh5OnlyKeepInvalid,featureDataColumnName= c("proteinName")) 

### Simply subtract channel 5 and 6, keep invalid
esetSubtractCh5and3KeepInvalid <-  eset
exprs(esetSubtractCh5and3KeepInvalid)[,5]  <- exprs(esetSubtractCh5and3KeepInvalid)[,5] - exprs(esetSubtractCh5and3KeepInvalid)[,4]*0.041 - exprs(esetSubtractCh5and3KeepInvalid)[,6]*0.033
exprs(esetSubtractCh5and3KeepInvalid)[,3]  <- exprs(esetSubtractCh5and3KeepInvalid)[,3] - exprs(esetSubtractCh5and3KeepInvalid)[,2]*0.05 - exprs(esetSubtractCh5and3KeepInvalid)[,4]*0.015
esetSubtractCh5and3KeepInvalid <- esetSubtractCh5and3KeepInvalid[apply(esetSubtractCh5and3KeepInvalid,1,function(t){ return(sum(is.na(t)) == 0 )}),]
esetProtSubtractCh5and3KeepInvalid <- rollUp(esetSubtractCh5and3KeepInvalid,featureDataColumnName= c("proteinName")) 

### standard purity correction, keep invalid
esetSolveCh5and3KeepInvalid <-  eset
exprs(esetSolveCh5and3KeepInvalid) <- purityCorrectTMT(exprs(esetSolveCh5and3KeepInvalid),getImpuritiesMatrix(-60), invalidReplace = "")
esetProtSolveCh5and3KeepInvalid <- rollUp(esetSolveCh5and3KeepInvalid,featureDataColumnName= c("proteinName")) 

### eset
#eset <-eset[,1:5]

#exprs(eset)[,5] <- exprs(eset)[,5] - exprs(eset)[,4]*0.041 - exprs(eset)[,6]*0.033


#naBefore <- sum(is.na(exprs(eset)) | exprs(eset) < 0  )
#exprs(eset)[,5] <- exprs(eset)[,5] - exprs(eset)[,4]*0.041 
#naAfter <- sum(is.na(exprs(eset)) | exprs(eset) < 0  )
#exprs(eset)[ exprs(eset) < 0] <- NA

#print(nrow(exprs(eset)) - sum(  apply(exprs(eset),1,function(t){sum(is.na(t) | (t < 0)) > 0 }) ))

#sel <- (exprs(eset)[,5] - exprs(eset)[,4]*0.041 - exprs(eset)[,6]*0.033) < 0
#exprs(eset)[exprs(eset) < 0 ] <- 1


### discard features with missing values
#eset <- eset[apply(eset,1,function(t){ return(sum(is.na(t)) == 0 )}),]


#pdf("/tmp/tmp.pdf")
#pdf("/tmp/tmp_pCorrect_diag1.pdf")
pdf("/tmp/tmp_pCorrect.pdf")
par(mfrow=c(2,2))
#pdf("/tmp/tmp_discard_interference_ch4_from_ch_5_manual.pdf")

selChannels <- 1:5

myRatioBoxplot(esetProtSolveAllZero[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"esetProtSolveAllZero")

myRatioBoxplot(esetProtSolveAllNA[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"esetProtSolveAllNA")

myRatioBoxplot(esetProtSolveAllOrg[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"esetProtSolveAllOrg")

myRatioBoxplot(esetProtSolveKeepInvalid[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"esetProtSolveKeepInvalid")

par(mfrow=c(2,2))
myRatioBoxplot(esetProtSubtractCh5OnlyKeepInvalid[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"SubtractCh5OnlyKeepInvalid")

myRatioBoxplot(esetProtSolveCh5OnlyKeepInvalid[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"SolveCh5OnlyKeepInvalid")

myRatioBoxplot(esetProtSubtractCh5and3KeepInvalid[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"SubtractCh5and3KeepInvalid")

myRatioBoxplot(esetProtSolveCh5and3KeepInvalid[,selChannels],dispN=T,cex=1)
mtext(side=3,at=1,"SolveCh5and3KeepInvalid")

dev.off()


#boxplot(exprs(eset),log="y")
