# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### load / source
library("affy")
library("limma")
library(gplots) # volcano plot
library(seqinr)
library(optparse)
library(data.table)

##@TEMP
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")

source("/Users/erikahrne/dev/R/workspace/TMTRatioCorrection/TMTRatioCorrection.R")
### PARAMS

#10 PLEX:
#RPE Channel 1&2
#Hela Channel 3&4
#DLD1 2N Channel 5
#DLD1 4N Channel 6
#DLD1 3N clone 7 Channel 7
#DLD1 3N clone 13 Channel 8
#DLD1 3N clone 16 Channel 9
#DLD1 more 4N clone 11 Channel 10

if(T){
	rep1ScaffoldFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20150326-091346_TMTex28b/Scaffold/Raw Data Report for ENigg-Cristine-TMT-28b-40-41-110515-human-short.xls"
	rep2ScaffoldFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20150326-091346_TMTex28b/Scaffold/Raw Data Report for ENigg-Cristina-TMT-28br-40r-41r-110515-human-short.xls"
	rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/caseStudyCristinaTriplicate10plexTMP_short.rData"

	#CTRL <- "RPE"
	#CTRL <- "DLD12N"
	CTRL <- "DLD14N"
	
	pdfFile <- paste("/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/short_",CTRL,".pdf",sep="")
	xlsFile <- paste("/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/short_",CTRL,".xls",sep="")

	
}else{
	rep1ScaffoldFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20150326-091346_TMTex28b/Scaffold/Raw Data Report for ENigg-Cristine-TMT-28b-40-41-no-replicates-200415-final.xls"
	rep2ScaffoldFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20150326-091346_TMTex28b/Scaffold/Raw Data Report for ENigg-Cristina-TMT-28br-40r-41r-no-replicates-200415-final.xls"
	pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/wIsoForms.pdf"
	xlsFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyCristinaTriplicate10Plex/wIsoForms.xls"
	rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/caseStudyCristinaTriplicate10plexTMP_wIsoForms.rData"
	
}

qValueThrs <- 0.05

### PARAMS END

condition <- c("RPE","RPE","HeLa","HeLa","DLD12N","DLD14N","DLD13Nclone7","DLD13Nclone13","DLD13Nclone16","DLD1more4Nclone11")
expDesign <- data.frame(condition=condition,isControl=grepl(CTRL,condition) )

if(F){
	print("PARSE")
	
	esetTechRep1 <- parseScaffoldRawFile(file=rep1ScaffoldFile,expDesign=expDesign)
	esetTechRep2 <- parseScaffoldRawFile(file=rep2ScaffoldFile,expDesign=expDesign)
	print("PARSE END")
	save(esetTechRep1,esetTechRep2,file=rDataFile)
}else{
	load(rDataFile)
}

# test @ TMP
#esetTechRep1 <- esetTechRep1[rownames(esetTechRep1) %in% sample(nrow(esetTechRep1),3000),] 
#esetTechRep2 <- esetTechRep2[rownames(esetTechRep2) %in% sample(nrow(esetTechRep2),3000),] 

### Merge technical replicate

expressionMatrix=rbind(exprs(esetTechRep1),exprs(esetTechRep2))
featureAnnotations=rbind(fData(esetTechRep1),fData(esetTechRep2))
rownames(expressionMatrix) <-1:nrow(expressionMatrix)
rownames(featureAnnotations) <-1:nrow(expressionMatrix)
eset <- createExpressionDataset(expressionMatrix,expDesign=expDesign,featureAnnotations)

# update expDesign
pData(eset)$isControl <- grepl(CTRL,pData(eset)$condition)

# discard Decoys and Contaminants
eset <- eset[!(isDecoy(as.character(fData(eset)$proteinName)) | isCon(as.character(fData(eset)$proteinName))),]


### Seperate biological replicates. Biological Replicate Labels 28b, 40, 41 
#	TMTex28b TMTex28b_150416091503             TMTex28br	TMTex28brr
#	53900                  4872                 25259		47721 		=  131752
#	TMTexp40             TMTexp40r 
#   55160                 79069 	= 134229
#	TMTexp41             TMTexp41r 
#	34620                 57787 = 92407

# biological replicate 1
#sum(grepl("ex28",as.character(fData(eset)$spectrumName)))
esetBioRep1 <- eset[grepl("ex28",as.character(fData(eset)$spectrumName)),]

# biological replicate 2
#sum(grepl("exp40",as.character(fData(eset)$spectrumName)))
esetBioRep2 <- eset[grepl("exp40",as.character(fData(eset)$spectrumName)),]

# biological replicate 3
#sum(grepl("exp41",as.character(fData(eset)$spectrumName)))
esetBioRep3 <- eset[grepl("exp41",as.character(fData(eset)$spectrumName)),]

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

# roll-up

# calculate ratios
proteinRatiosBioRep1 <- getRatios(esetBioRep1ProteinNorm)
proteinRatiosBioRep2 <- getRatios(esetBioRep2ProteinNorm)
proteinRatiosBioRep3 <- getRatios(esetBioRep3ProteinNorm)

#tmp <- esetBioRep1ProteinNorm[1:10,]
#pData(tmp)$isControl <- c(F,F,F,F,T,F,F,F,F,F)
#getRatios(tmp)

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


expDesign <- data.frame(row.names=1:21,isControl=rep(F,21),condition=unlist(lapply(names(proteinRatiosBioRep1),function(t){ rep(t,3)})))
#expDesign <- data.frame(row.names=1:21,isControl=rep(F,21),condition=   paste(1:3,names(proteinRatiosBioRep1))   )
esetProteinRatios <- createExpressionDataset(expressionMatrix=ratioExpMatrix,expDesign=expDesign,featureAnnotations=featureData)		 

###### MERGE DATA END

###### DIFFERENTIAL ABUNDANCE STATS
pValues <- data.frame(eBayes(lmFit(esetProteinRatios[,1:3]))$p.value
		,eBayes(lmFit(esetProteinRatios[,4:6]))$p.value
		,eBayes(lmFit(esetProteinRatios[,7:9]))$p.value
		,eBayes(lmFit(esetProteinRatios[,10:12]))$p.value
		,eBayes(lmFit(esetProteinRatios[,13:15]))$p.value
		,eBayes(lmFit(esetProteinRatios[,16:18]))$p.value
		,eBayes(lmFit(esetProteinRatios[,19:21]))$p.value
)
#names(pValues) <- paste("condition",2:8,sep="_")
names(pValues) <- names(proteinRatiosBioRep1)
qValues <- data.frame(apply(pValues,2,p.adjust,method="BH"))		
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

medianRatiosCalibrated <- data.frame(predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,1]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,2]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,3]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,4]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,5]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,6]))
									,predict(fitRatioCalibrationProtein,newdata=data.frame(tmtRatio=medianRatios[,7]))
									,row.names=rownames(medianRatios)		
							)

names(medianRatiosCalibrated)	<- names(medianRatios)						
							
### ADJUST RATIOS BASED ON SPIKE-IN KIT END

### XLS

ratios <- exprs(esetProteinRatios)
colnames(ratios) <- paste(pData(esetProteinRatios)$condition ,colnames(exprs(esetProteinRatios)),sep="_")
out <- cbind(fData(esetProteinRatios),ratios,medianRatios,medianRatiosCalibrated,qValues)
names(out)[3:23] <- paste("log2ratio",names(out)[3:23],sep="_")
names(out)[24:30] <- paste("median_log2ratio",names(out)[24:30],sep="_")
names(out)[31:37] <- paste("median_log2ratio_adjusted",names(out)[31:37],sep="_")
names(out)[38:ncol(out)] <- paste("qValue",names(out)[38:ncol(out)],sep="_")
write.table(file=xlsFile,out,row.names=F,sep="\t")

### XLS END

### GRAPHICS

pdf(pdfFile)

par(mfrow=c(2,3))
plotCalibrationMixRatios(esetCalibMix)
par(mfrow=c(1,1))
plotRelativeInterferencePerChannel(esetCalibMix, cex.axis=1.5, cex.lab=1.5)
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

par(mfrow=c(1,2))
venn(list(bioRep1=unique(as.character(fData(esetBioRep1)$peptide)),bioRep2=unique(as.character(fData(esetBioRep2)$peptide)),bioRep3=unique(as.character(fData(esetBioRep3)$peptide))))
mtext(side=3,"Peptides")

venn(list(bioRep1=unique(as.character(fData(esetBioRep1)$proteinName)),bioRep2=unique(as.character(fData(esetBioRep2)$proteinName)),bioRep3=unique(as.character(fData(esetBioRep3)$proteinName))))
mtext(side=3,"Proteins")

par(mfrow=c(3,3),cex.lab=1.5,cex.names=1.5, cex.axis=1.3 )
for(i in 1:ncol(qValues)){
	#controlCondition <- "RPE"
	#caseCondition <- paste("cond", i+1)
	caseCondition <- names(medianRatios)[i]		
	
	plot(medianRatios[,i],-log10(qValues[,i])
			, pch=19
			, xlab= paste("log2(",caseCondition,"/",CTRL,")", sep="" )
			, ylab= paste("-log10(qValue)",sep="")
			,col=c("black","blue")[(qValues[,i] < qValueThrs)+1 ]
	
	
	)
	grid()
}

barplot(apply(qValues <= qValueThrs,2,sum  ),las=2, names=names(medianRatios), ylab=paste("Diff Abd. ( qValue <=",qValueThrs,")"))


dev.off()

### GRAPHICS END


print("DONE")