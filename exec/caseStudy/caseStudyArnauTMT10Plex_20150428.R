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

scaffoldRawFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MZavolan/ArnauVinaVilaseca_14/20150428-172535_p40/Scaffold/Raw Data Report for MZavolan-Arnau-TMT-p40-040615.xls"



rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/caseStudyArnauTMT10Plex_20150428.rData"

CTRL <- "l4"
qValueThrs <- 0.01

pdfFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MZavolan/ArnauVinaVilaseca_14/20150428-172535_p40/Scaffold/SafeQuant/SQ.pdf"
xlsFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MZavolan/ArnauVinaVilaseca_14/20150428-172535_p40/Scaffold/SafeQuant/SQ.xls"

### PARAMS END

### PARSE

# 	"passage rep1" "passage rep2" "H/L" "condition rep1"
#126	4	4	H	h4	h4
#127N	4	4	L	l4	l4
#127C	5	4	H	h5	h4
#128N	5	4	L	l5	l4
#128C	4	5	H	h4	h5
#129N	4	5	L	l4 	l5
#129C	5	5	H	h5	h5
#130N	5	5	L	l5	l5
#130C	6	6	H	h6	h6
#131	6	6	L	l6	l6

#conditions high4, low4, high5, low5, high6, low6 

if(F){
	print("PARSE")
	
	expDesign <- data.frame(condition=paste("Condition",c(1,2,3,1,2,3,1,2,3,1)),isControl=c(T,F,F,T,F,F,T,F,F,T) )
	esetTmp <- parseScaffoldRawFile(file=scaffoldRawFile,expDesign=expDesign)
	esetTmp <- esetTmp[!isDecoy(fData(esetTmp)$proteinName) & !isCon(fData(esetTmp)$proteinName),]
	print("PARSE END")
	save(esetTmp,file=rDataFile)
}else{
	load(rDataFile)
}

fData(esetTmp)$spectrumName <- as.character(fData(esetTmp)$spectrumName)
runNb <- as.numeric(gsub("_p-40.*","",fData(esetTmp)$spectrumName))

eset1 <- esetTmp[runNb < 11,]
eset2 <- esetTmp[runNb > 10,]

pData(eset1)$condition <- c("h4_d1_1","l4_d1_1","h5_d1_1","l5_d1_1","h4_d2_1","l4_d2_1","h5_d2_1","l5_d2_1","h6_d2_1","l6_d2_1")
pData(eset1)$isControl <- grepl(CTRL,pData(eset1)$condition)
pData(eset1)$condition[pData(eset1)$isControl] <- gsub("\\_.*","",pData(eset1)$condition[pData(eset1)$isControl])
pData(eset1)$condition <- as.factor(pData(eset1)$condition)


pData(eset2)$condition <- c("h4_d1_2","l4_d1_2","h4_d3_2","l4_d3_2","h5_d1_2","l5_d1_2","h5_d3_2","l5_d3_2","h6_d3_2","l6_d3_2")
pData(eset2)$isControl <- grepl(CTRL,pData(eset2)$condition) 
pData(eset2)$condition[pData(eset2)$isControl] <- gsub("\\_.*","",pData(eset2)$condition[pData(eset2)$isControl])
pData(eset2)$condition <- as.factor(pData(eset2)$condition)

pData(eset2)



### PARSE END

### PRE-PROCESS

# normalize
# discard runs
eset1Norm <- eset1[,2:10]
eset2Norm <- eset2[,2:10]

eset1Norm <- sqNormalize(eset1, method="global")
eset2Norm <- sqNormalize(eset2, method="global")

# replace missing values
baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(esetTmp)[,1])),promille=5)

exprs(eset1Norm)[is.na(exprs(eset1Norm)) | (exprs(eset1Norm) < 0)  ] <- 0 
exprs(eset1Norm) <- exprs(eset1Norm) + baselineIntensity

pData(eset1Norm)
# pData(eset1Norm)
#    condition isControl globalNormFactors
# 1    h4_d1_1     FALSE        1.00000000
# 2         l4      TRUE        0.13008807
# 3    h5_d1_1     FALSE        0.10134815
# 4    l5_d1_1     FALSE        0.11227925
# 5    h4_d2_1     FALSE        0.06774442
# 6         l4      TRUE        0.09846373
# 7    h5_d2_1     FALSE        0.09836274
# 8    l5_d2_1     FALSE        0.10973232
# 9    h6_d2_1     FALSE        0.10909501
# 10   l6_d2_1     FALSE        0.10640477
# > cat("Synch1435221768638116000\n");


exprs(eset2Norm)[is.na(exprs(eset2Norm)) | (exprs(eset2Norm) < 0)  ] <- 0 
exprs(eset2Norm) <- exprs(eset2Norm) + baselineIntensity

pData(eset2Norm)
# pData(eset2Norm)
#    condition isControl globalNormFactors
# 1    h4_d1_2     FALSE        1.00000000
# 2         l4      TRUE        0.12657107
# 3    h4_d3_2     FALSE        0.15312200
# 4         l4      TRUE        0.10409090
# 5    h5_d1_2     FALSE        0.06645898
# 6    l5_d1_2     FALSE        0.09392171
# 7    h5_d3_2     FALSE        0.13872919
# 8    l5_d3_2     FALSE        0.11540657
# 9    h6_d3_2     FALSE        0.10374254
# 10   l6_d3_2     FALSE        0.08047145
# > cat("Synch1435221908138571000\n");



print("rollUp")
eset1ProteinNorm <- rollUp(eset1Norm)
eset2ProteinNorm <- rollUp(eset2Norm)
print("rollUp DONE")


### PRE-PROCESS END

# calculate ratios
proteinRatios1 <- getRatios(eset1ProteinNorm)
proteinRatios2 <- getRatios(eset2ProteinNorm)

# names(proteinRatios1)
# [1] "h4_d1_1" "h5_d1_1" "l5_d1_1" "h4_d2_1" "h5_d2_1" "l5_d2_1" "h6_d2_1"
# [8] "l6_d2_1"

# names(proteinRatios2)
# [1] "h4_d1_2" "h4_d3_2" "h5_d1_2" "l5_d1_2" "h5_d3_2" "l5_d3_2" "h6_d3_2"
# [8] "l6_d3_2"


###### MERGE DATA

allProteins <-  unique(c(rownames(eset1ProteinNorm),rownames(eset2ProteinNorm)))

ratioExpMatrix <- cbind(proteinRatios1[allProteins,],proteinRatios2[allProteins,])
ratioExpMatrix <- ratioExpMatrix[,sort(names(ratioExpMatrix))]
rownames(ratioExpMatrix) <- allProteins

featureData <- unique(rbind(fData(eset1ProteinNorm)[,1:2]
				,fData(eset2ProteinNorm)[,1:2]
				
		))[allProteins,]



expDesign <- data.frame(row.names=names(ratioExpMatrix),isControl=rep(F,ncol(ratioExpMatrix)),condition=gsub("\\_.*","",names(ratioExpMatrix)) )
esetProteinRatios <- createExpressionDataset(expressionMatrix=as.matrix(ratioExpMatrix),expDesign=expDesign,featureAnnotations=featureData)		 

pData(esetProteinRatios)
# pData(esetProteinRatios)
#         isControl condition
# h4_d1_1     FALSE        h4
# h4_d1_2     FALSE        h4
# h4_d2_1     FALSE        h4
# h4_d3_2     FALSE        h4
# h5_d1_1     FALSE        h5
# h5_d1_2     FALSE        h5
# h5_d2_1     FALSE        h5
# h5_d3_2     FALSE        h5
# h6_d2_1     FALSE        h6
# h6_d3_2     FALSE        h6
# l5_d1_1     FALSE        l5
# l5_d1_2     FALSE        l5
# l5_d2_1     FALSE        l5
# l5_d3_2     FALSE        l5
# l6_d2_1     FALSE        l6
# l6_d3_2     FALSE        l6
# > cat("Synch1435225504006630000\n");


###### MERGE DATA END

###### DIFFERENTIAL ABUNDANCE STATS

# pData(esetProteinRatios)

pValues <- data.frame(eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition == levels(pData(esetProteinRatios)$condition)[1] ]))$p.value
	,eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition ==  levels(pData(esetProteinRatios)$condition)[2] ]))$p.value
	,eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition == levels(pData(esetProteinRatios)$condition)[3] ]))$p.value
	,eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition == levels(pData(esetProteinRatios)$condition)[4] ]))$p.value
	,eBayes(lmFit(esetProteinRatios[,pData(esetProteinRatios)$condition == levels(pData(esetProteinRatios)$condition)[5] ]))$p.value
	
)
#names(pValues) <- paste("condition",2:8,sep="_")
medianRatios <- getSignalPerCondition(esetProteinRatios)

names(pValues) <- names(medianRatios)
qValues <- data.frame(apply(pValues,2,p.adjust,method="BH"))		

sqa <- list()
class(sqa) <- "safeQuantAnalysis"

sqa$eset <- esetProteinRatios
sqa$qValue <- qValues
sqa$pValue <- pValues
sqa$ratio <- medianRatios

###### DIFFERENTIAL ABUNDANCE STATS END


### GRAPHICS

pdf(pdfFile)

plotExpDesign(eset1)
plotExpDesign(eset2)

par(mfrow=c(1,2))
barplotMSSignal(eset1, main="set 1", labels=pData(eset1)$condition)
barplotMSSignal(eset2, main="set 2", labels=pData(eset2)$condition)

par(mfrow=c(1,1))

conditionColors =.getConditionColors(esetProteinRatios)
samplecolors =  as.vector(unlist(conditionColors[pData(esetProteinRatios)$condition,]))
breaks=seq(-1.5,1.5,length=20)

m <- exprs(esetProteinRatios)
m[is.na(m)] <- 0

heatmap.2(
		m
		,ColSideColors=samplecolors
		,breaks=breaks
		, col=colorRampPalette(c(colors()[142],"black",colors()[128]))
		,labRow= rep("",(nrow(m)))
		,dendrogram="column"
		,density.info="density", scale="none"
		, key=TRUE
		, symkey=FALSE
		, trace="none"
)
legend("left",levels(pData(esetProteinRatios)$condition), fill=as.character(conditionColors[,1]), cex=0.7, box.col=0)

if(F){
	# volcanoes qvalue
	par(mfrow=c(3,3),cex.lab=1.5,cex.axis=1.3 )
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
	par(mfrow=c(3,3),cex.lab=1.5, cex.axis=1.3 )
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
}

#barplot(apply(qValues <= qValueThrs,2,sum  ),las=2, names=names(medianRatios), ylab=paste("Diff Abd. ( qValue <=",qValueThrs,")"))
#barplot(apply(pValues <= qValueThrs,2,sum  ),las=2, names=names(medianRatios), ylab=paste("Diff Abd. ( pValue <=",qValueThrs,")"))

par(mfrow=c(2,2))
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=T, isLegend=T ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=T, isLegend=F, pvalCutOff=qValueThrs)

plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=F , isLegend=F ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=F, isLegend=F ,pvalCutOff=qValueThrs)

par(mfrow=c(3,3))
.allpValueHist(sqa)

dev.off()

### GRAPHICS END


### XLS

ratios <- exprs(esetProteinRatios)
colnames(ratios) <- paste("log2ratio",colnames(exprs(esetProteinRatios)),sep="_")
out <- cbind(fData(esetProteinRatios),ratios,medianRatios,pValues,qValues)
names(out)[19:23] <- paste("median_log2ratio",names(out)[19:23],sep="_")

names(out)[24:28] <- paste("pValue",names(out)[24:28],sep="_")
names(out)[29:ncol(out)] <- paste("qValue",names(out)[29:ncol(out)],sep="_")
write.table(file=xlsFile,out,row.names=F,sep="\t")

### XLS END


print("DONE")