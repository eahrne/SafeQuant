# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### load / source
library("affy")
library("limma")
library(gplots) # volcano plot
library(seqinr)
library(optparse)
library(data.table)

##@TEMP
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/TMT.R")
source("/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/UserOptions.R")

source("/Users/ahrnee-adm/dev/R/workspace/TMTRatioCorrection/TMTRatioCorrection.R")
### PARAMS


scaffoldRawFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20151005-173040_exp64/Scaffold/Raw Data Report for Copy_of_ENigg-Christina-TMT-Merge-exp62-63-64-28102015.xls" 
rDataFile <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyTMTMergeBioReplicates/caseStudyTMTMergeBioReplicates.rData"
#pdfFile <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/20151112/exp_62-64/exp_62-64.pdf"
#xlsFile <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/20151112/exp_62-64/exp_62-64.xls"

pdfFile <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/temp/exp_62-64.pdf"
xlsFile <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/temp/exp_62-64.xls"

#outSel <- 1:21 # 
#outSel <- 13:21 
outSel <- 1:9

### PARAMS END

expDesign <- data.frame(condition=paste("Condition",c(1,1,2,3,4,5,5,6,7,8)),isControl=c(T,T,F,F,F,F,F,F,F,F) )
#expDesign <- data.frame(condition=paste("Condition",c(1,1,2,3,4,5,5,6,7,8)),isControl=c(F,F,F,F,F,T,T,F,F,F) )
if(T){
	print("PARSE")
	
	eset <- parseScaffoldRawFile(file=scaffoldRawFile,expDesign=expDesign)
	print("PARSE END")
	save(eset,file=rDataFile)
}else{
	load(rDataFile)
}

pData(eset) <-  expDesign
CTRL <- expDesign$condition[expDesign$isControl][1]
qValueThrs <- 0.05

### add CTRL tag to output files
pdfFile <- gsub(" ","_",gsub("\\.pdf",paste("_",CTRL,".pdf",sep=""),pdfFile))
xlsFile <- gsub(" ","_",gsub("\\.xls",paste("_",CTRL,".xls",sep=""),xlsFile))

# discard Decoys and Contaminants
eset <- eset[!(isDecoy(as.character(fData(eset)$proteinName)) | isCon(as.character(fData(eset)$proteinName))),]

# biological replicate 1
esetBioRep1 <- eset[grepl("exp\\-62",as.character(fData(eset)$spectrumName)),]

# biological replicate 2
esetBioRep2 <- eset[grepl("exp\\-63",as.character(fData(eset)$spectrumName)),]

# biological replicate 3
esetBioRep3 <- eset[grepl("exp\\-64",as.character(fData(eset)$spectrumName)),]

# keep only proteins found in all biological replicates
sharedProteins <- unique(intersect(intersect(fData(esetBioRep1)$proteinName,fData(esetBioRep2)$proteinName),fData(esetBioRep3)$proteinName))
esetBioRep1 <- esetBioRep1[fData(esetBioRep1)$proteinName %in% sharedProteins, ]
esetBioRep2 <- esetBioRep2[fData(esetBioRep2)$proteinName %in% sharedProteins, ]
esetBioRep3 <- esetBioRep3[fData(esetBioRep3)$proteinName %in% sharedProteins, ]

# normalize
esetBioRep1Norm <- sqNormalize(esetBioRep1, method="global")
esetBioRep2Norm <- sqNormalize(esetBioRep2, method="global")
esetBioRep3Norm <- sqNormalize(esetBioRep3, method="global")

# replace missing values
baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(eset)[,1])),promille=5)

exprs(esetBioRep1Norm)[is.na(exprs(esetBioRep1Norm)) | (exprs(esetBioRep1Norm) < 0)  ] <- 0 
exprs(esetBioRep1Norm) <- exprs(esetBioRep1Norm) + baselineIntensity

exprs(esetBioRep2Norm)[is.na(exprs(esetBioRep2Norm)) | (exprs(esetBioRep2Norm) < 0)  ] <- 0 
exprs(esetBioRep2Norm) <- exprs(esetBioRep2Norm) + baselineIntensity

exprs(esetBioRep3Norm)[is.na(exprs(esetBioRep3Norm)) | (exprs(esetBioRep3Norm) < 0)  ] <- 0 
exprs(esetBioRep3Norm) <- exprs(esetBioRep3Norm) + baselineIntensity

print("rollUp")
esetBioRep1ProteinNorm <- rollUp(esetBioRep1Norm)
esetBioRep2ProteinNorm <- rollUp(esetBioRep2Norm)
esetBioRep3ProteinNorm <- rollUp(esetBioRep3Norm)

print("rollUp DONE")

# calculate ratios
proteinRatiosBioRep1 <- getRatios(esetBioRep1ProteinNorm)
proteinRatiosBioRep2 <- getRatios(esetBioRep2ProteinNorm)
proteinRatiosBioRep3 <- getRatios(esetBioRep3ProteinNorm)

###### MERGE DATA

allProteins <-  unique(c(rownames(proteinRatiosBioRep1),rownames(proteinRatiosBioRep2),rownames(proteinRatiosBioRep3)))
ratioExpMatrix <- cbind(proteinRatiosBioRep1[allProteins,1], proteinRatiosBioRep2[allProteins,1] ,proteinRatiosBioRep3[allProteins,1]
		,proteinRatiosBioRep1[allProteins,2], proteinRatiosBioRep2[allProteins,2] ,proteinRatiosBioRep3[allProteins,2]
		,proteinRatiosBioRep1[allProteins,3], proteinRatiosBioRep2[allProteins,3] ,proteinRatiosBioRep3[allProteins,3]
		,proteinRatiosBioRep1[allProteins,4], proteinRatiosBioRep2[allProteins,4] ,proteinRatiosBioRep3[allProteins,4]
		,proteinRatiosBioRep1[allProteins,5], proteinRatiosBioRep2[allProteins,5] ,proteinRatiosBioRep3[allProteins,5]
		,proteinRatiosBioRep1[allProteins,6], proteinRatiosBioRep2[allProteins,6] ,proteinRatiosBioRep3[allProteins,6]
		,proteinRatiosBioRep1[allProteins,7], proteinRatiosBioRep2[allProteins,7] ,proteinRatiosBioRep3[allProteins,7]
)

featureData <- unique(rbind(fData(esetBioRep1ProteinNorm)[,1:2]
				,fData(esetBioRep2ProteinNorm)[,1:2]
				,fData(esetBioRep3ProteinNorm)[,1:2]
		))[allProteins,]

expDesignMerge <- data.frame(row.names=1:ncol(ratioExpMatrix),isControl=rep(F,ncol(ratioExpMatrix)),condition=unlist(lapply(names(proteinRatiosBioRep1),function(t){ rep(t,3)})))
#expDesign <- data.frame(row.names=1:21,isControl=rep(F,21),condition=   paste(1:3,names(proteinRatiosBioRep1))   )
esetProteinRatios <- createExpressionDataset(expressionMatrix=ratioExpMatrix,expDesign=expDesignMerge,featureAnnotations=featureData)		 

### filter
esetProteinRatios <- esetProteinRatios[,outSel]


###### MERGE DATA END

###### DIFFERENTIAL ABUNDANCE STATS
#pValues <- data.frame(eBayes(lmFit(esetProteinRatios[,1:3]))
#		,eBayes(lmFit(esetProteinRatios[,4:6]))$p.value
#		,eBayes(lmFit(esetProteinRatios[,7:9]))$p.value
#		,eBayes(lmFit(esetProteinRatios[,10:12]))$p.value
#		,eBayes(lmFit(esetProteinRatios[,13:15]))$p.value
#		,eBayes(lmFit(esetProteinRatios[,16:18]))$p.value
#		,eBayes(lmFit(esetProteinRatios[,19:21]))$p.value
#)
pValues <- data.frame(row.names=allProteins)
for(cond in unique(pData(esetProteinRatios)$condition)){
	pValues <- cbind(pValues, eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition == cond]))$p.value )
}

#names(pValues) <- paste("condition",2:8,sep="_")
names(pValues) <- unique(pData(esetProteinRatios)$condition)
qValues <- data.frame(apply(pValues,2,p.adjust,method="BH"),check.names=F)	
medianRatios <- getSignalPerCondition(esetProteinRatios)
###### DIFFERENTIAL ABUNDANCE STATS END

### ADJUST RATIOS BASED ON SPIKE-IN KIT

esetCalibMix <- .getCalibMixEset(eset)
esetCalibMixPair <- .getCalibMixPairedEset(esetCalibMix)
fitRatioCalibrationProtein <- getRatioCorrectionFactorModel(rollUp(esetCalibMixPair,featureDataColumnName="proteinName"))
#fitRatioCalibrationPeptide <- getRatioCorrectionFactorModel(rollUp(esetCalibMixPair,featureDataColumnName="peptide"))
#fitRatioCalibrationSpectrum <- getRatioCorrectionFactorModel(esetCalibMixPair)

#getRatioCorrectionFactorModel(rollUp(.getCalibMixPairedEset(.getCalibMixEset(eset)),featureDataColumnName="peptide"))
#getRatioCorrectionFactorModel(rollUp(.getCalibMixPairedEset(.getCalibMixEset(eset)),featureDataColumnName="proteinName"))
#getRatioCorrectionFactorModel(rollUp(.getCalibMixPairedEset(.getCalibMixEset(esetBioRep1)),featureDataColumnName="proteinName"))
#getRatioCorrectionFactorModel(rollUp(.getCalibMixPairedEset(.getCalibMixEset(esetBioRep2)),featureDataColumnName="proteinName"))
#getRatioCorrectionFactorModel(rollUp(.getCalibMixPairedEset(.getCalibMixEset(esetBioRep3)),featureDataColumnName="proteinName"))

#medianRatiosCalibrated <- data.frame(predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,1]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,2]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,3]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,4]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,5]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,6]))
#									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,7]))
#									,row.names=rownames(medianRatios)		
#							)

print(fitRatioCalibrationProtein)

medianRatiosCalibrated <-data.frame(row.names=allProteins)
for(i in 1:ncol(medianRatios)){
	medianRatiosCalibrated <- cbind(medianRatiosCalibrated,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,i])))
}							
names(medianRatiosCalibrated)	<- names(medianRatios)						

sqa <- list()
class(sqa) <- "safeQuantAnalysis"

sqa$eset <- esetProteinRatios
sqa$qValue <- qValues
sqa$pValue <- pValues
sqa$ratio <- medianRatios

### ADJUST RATIOS BASED ON SPIKE-IN KIT END

### GRAPHICS

pdf(pdfFile)

conditionColors =.getConditionColors(esetProteinRatios)
samplecolors =  as.vector(unlist(conditionColors[pData(esetProteinRatios)$condition,]))
breaks=seq(-1.5,1.5,length=20)

heatmap.2(exprs(esetProteinRatios)
		
		,ColSideColors=samplecolors
		,breaks=breaks
		, col=colorRampPalette(c(colors()[142],"black",colors()[128]))
		,labRow= rep("",(nrow(esetProteinRatios)))
		,dendrogram="column"
		,density.info="density", scale="none"
		, key=TRUE
		, symkey=FALSE
		, trace="none"
)
legend("left",levels(as.factor(pData(esetProteinRatios)$condition)), fill=as.character(conditionColors[,1]), cex=0.7, box.col=0)

#par(mfrow=c(2,3))
#plotCalibrationMixRatios(esetCalibMix)
#par(mfrow=c(1,1))
# plotRelativeInterferencePerChannel(esetCalibMix, cex.axis=1.5, cex.lab=1.5)
#.plotTMTRatioVsRefRatio(esetCalibMixPair)
#.plotTMTRatioVsRefRatio(rollUp(esetCalibMixPair,featureDataColumnName="proteinName"))
.plotTMTRatioVsRefRatio(rollUp(esetCalibMixPair,featureDataColumnName="peptide"))

par(mfrow=c(2,2))
missinValueBarplot(esetBioRep1, main="BioRep 1")
missinValueBarplot(esetBioRep2, main="BioRep 2")
missinValueBarplot(esetBioRep3, main="BioRep 3")

par(mfrow=c(2,2))
barplotMSSignal(esetBioRep1, main="BioRep 1")
barplotMSSignal(esetBioRep2, main="BioRep 2")
barplotMSSignal(esetBioRep3, main="BioRep 3")

#par(mfrow=c(1,2))
#venn(list(bioRep1=unique(as.character(fData(esetBioRep1)$peptide)),bioRep2=unique(as.character(fData(esetBioRep2)$peptide)),bioRep3=unique(as.character(fData(esetBioRep3)$peptide))))
#mtext(side=3,"Peptides")
#
#venn(list(bioRep1=unique(as.character(fData(esetBioRep1)$proteinName)),bioRep2=unique(as.character(fData(esetBioRep2)$proteinName)),bioRep3=unique(as.character(fData(esetBioRep3)$proteinName))))
#mtext(side=3,"Proteins")

# qvalue
par(mfrow=c(3,3),cex.lab=1.5,cex.names=1.5, cex.axis=1.3 )
for(i in 1:ncol(qValues)){
	caseCondition <- names(medianRatios)[i]		
	plot(medianRatios[,i],-log10(qValues[,i])
			, pch=19
			, xlab= paste("log2(",caseCondition,"/",CTRL,")", sep="" )
			, ylab= paste("-log10(qValue)",sep="")
			,col=c("black","blue")[(qValues[,i] < qValueThrs)+1 ]
	)
	grid()
}

# pvalue
par(mfrow=c(3,3),cex.lab=1.5,cex.names=1.5, cex.axis=1.3 )
for(i in 1:ncol(qValues)){
	caseCondition <- names(medianRatios)[i]		
	plot(medianRatios[,i],-log10(pValues[,i])
			, pch=19
			, xlab= paste("log2(",caseCondition,"/",CTRL,")", sep="" )
			, ylab= paste("-log10(pValue)",sep="")
			,col=c("black","blue")[(pValues[,i] < qValueThrs)+1 ]
	)
	grid()
}

#par(mfrow=c(2,2))
#barplot(apply(qValues <= qValueThrs,2,sum  ),las=2, names=names(medianRatios), ylab=paste("Diff Abd. ( qValue <=",qValueThrs,")"))
#barplot(apply(pValues <= qValueThrs,2,sum  ),las=2, names=names(medianRatios), ylab=paste("Diff Abd. ( pValue <=",qValueThrs,")"))

par(mfrow=c(2,2))
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=T, isLegend=T ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=T, isLegend=F, pvalCutOff=qValueThrs)

#pdf("/tmp/tmp.pdf")
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=F , isLegend=F ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=F, isLegend=F ,pvalCutOff=qValueThrs)
#dev.off()

#par(mfrow=c(3,3))
#.allpValueHist(sqa)

dev.off()

### GRAPHICS END

### XLS

ratios <- exprs(esetProteinRatios)
colnames(ratios) <- paste("log2ratio",paste(pData(esetProteinRatios)$condition ,colnames(exprs(esetProteinRatios)),sep="_"),sep="_")
colnames(medianRatios) <- paste("median_log2ratio",names(medianRatios),sep="_")
colnames(medianRatiosCalibrated) <- paste("median_log2ratio_adjusted",names(medianRatiosCalibrated),sep="_")
colnames(pValues) <- paste("pValue",names(pValues),sep="_")
colnames(qValues) <- paste("qValue",names(qValues),sep="_")
out <- cbind(fData(esetProteinRatios),ratios,medianRatios,medianRatiosCalibrated,pValues,qValues)
#names(out)[3:23] <- paste("log2ratio",names(out)[3:23],sep="_")
#names(out)[24:30] <- paste("median_log2ratio",names(out)[24:30],sep="_")
#names(out)[31:37] <- paste("median_log2ratio_adjusted",names(out)[31:37],sep="_")
#names(out)[38:44] <- paste("pValue",names(out)[38:44],sep="_")
#names(out)[45:ncol(out)] <- paste("qValue",names(out)[45:ncol(out)],sep="_")
write.table(file=xlsFile,out,row.names=F,sep="\t")

plot(medianRatios[,1],medianRatiosCalibrated[,1])
abline(coef=c(0,1))

### XLS END

print("DONE")