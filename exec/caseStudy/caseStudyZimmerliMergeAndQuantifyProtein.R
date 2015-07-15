# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### load / source
library("affy")
library("limma")
library(gplots)
library(seqinr)
library(optparse)
library(data.table)

##@TEMP
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")
#source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")
library(SafeQuant)
### PARAMS

progenesisFileRep1 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Proteome_Wide/Proteome_Wide_Set1.csv"
progenesisFileRep2 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Proteome_Wide/Proteome_Wide_Set2.csv"
progenesisFileRep3 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Proteome_Wide/Proteome_Wide_Set3.csv"
progenesisFileRep4 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Proteome_Wide/Proteome_Wide_Set4.csv"

rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/caseStudyZimmerliMergeAndQunatifyProtein.rData"

#pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyProtein/sq_intersect.pdf"
#xlsFile <-  "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyProtein/sq_intersect.xls"

pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyProtein/sq.pdf"
xlsFile <-  "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyProtein/sq.xls"

qValueThrs <- 0.01

### PARAMS END 
## parse
if(T){

	### PARSE
	
	esetR1 <- parseProgenesisPeptideMeasurementCsv(progenesisFileRep1,getExpDesignProgenesisCsv(progenesisFileRep1))
	esetR2 <- parseProgenesisPeptideMeasurementCsv(progenesisFileRep2,getExpDesignProgenesisCsv(progenesisFileRep2))
	esetR3 <- parseProgenesisPeptideMeasurementCsv(progenesisFileRep3,getExpDesignProgenesisCsv(progenesisFileRep3))
	esetR4 <- parseProgenesisPeptideMeasurementCsv(progenesisFileRep4,getExpDesignProgenesisCsv(progenesisFileRep4))
	
	### PARSE END
	
	### ROLL-UP
	
	esetProteinR1 <- rollUp(esetR1)
	esetProteinR2 <- rollUp(esetR2)
	esetProteinR3 <- rollUp(esetR3)
	esetProteinR4 <- rollUp(esetR4)
	
	### ROLL-UP END
	
	### SQA
	
	sqaProtR1 <- safeQuantAnalysis(esetProteinR1)
	sqaProtR2 <- safeQuantAnalysis(esetProteinR2)
	sqaProtR3 <- safeQuantAnalysis(esetProteinR3)
	sqaProtR4 <- safeQuantAnalysis(esetProteinR4)
	
	### SQA END 
		
	save(esetR1,esetR2,esetR3,esetR4,esetProteinR1,esetProteinR2,esetProteinR3,esetProteinR4,sqaProtR1,sqaProtR2,sqaProtR3,sqaProtR4,file=rDataFile)
	
}else{
	
	load(rDataFile)
	
}


allProteins <- unique(c( as.character(fData(sqaProtR1$eset)$proteinName)
						,as.character(fData(sqaProtR2$eset)$proteinName)
						,as.character(fData(sqaProtR3$eset)$proteinName)
						,as.character(fData(sqaProtR4$eset)$proteinName))
					)

		
intersectProteins <- intersect(
						intersect(as.character(fData(sqaProtR1$eset)$proteinName),as.character(fData(sqaProtR2$eset)$proteinName))
						,intersect(as.character(fData(sqaProtR3$eset)$proteinName),as.character(fData(sqaProtR4$eset)$proteinName))
				)
				
#selProteins	<- 	intersectProteins		
selProteins	<- 	allProteins		

expMatrixRatios <- data.frame(sqaProtR1$ratio[selProteins,],sqaProtR2$ratio[selProteins,],sqaProtR3$ratio[selProteins,],sqaProtR4$ratio[selProteins,])
rownames(expMatrixRatios) <- selProteins

expMatrixRatios <- expMatrixRatios[,rev(order(names(expMatrixRatios)))]

featureData <- unique(rbind(fData(sqaProtR1$eset)[,1:2]
				,fData(sqaProtR2$eset)[,1:2]
				,fData(sqaProtR3$eset)[,1:2]
				,fData(sqaProtR4$eset)[,1:2]
		
		))[selProteins,]

expDesign <- data.frame(row.names=names(expMatrixRatios),isControl=rep(F,ncol(expMatrixRatios)),condition=gsub("\\_[0-9]{1,}","",names(expMatrixRatios)) )
esetProteinRatios <- createExpressionDataset(expressionMatrix=as.matrix(expMatrixRatios),expDesign=expDesign,featureAnnotations=featureData)		 

#pData(esetProteinRatios)
#isControl condition
#WT_Rapa_155      FALSE   WT_Rapa
#WT_Rapa_1022     FALSE   WT_Rapa
#WT_Rapa_1020     FALSE   WT_Rapa
#WT_Rapa          FALSE   WT_Rapa
#KO_Rapa_150      FALSE   KO_Rapa
#KO_Rapa_1021     FALSE   KO_Rapa
#KO_Rapa_1012     FALSE   KO_Rapa
#KO_Rapa          FALSE   KO_Rapa
#KO_980           FALSE        KO
#KO_979           FALSE        KO
#KO_1104          FALSE        KO
#KO               FALSE        KO

# caluc median ratios per condition (merged data)
medianRatios <- getSignalPerCondition(esetProteinRatios)

# p-values
pValues <- data.frame(eBayes(lmFit(esetProteinRatios[,1:4 ]))$p.value
		,eBayes(lmFit(esetProteinRatios[,5:8 ]))$p.value
		,eBayes(lmFit(esetProteinRatios[,9:12 ]))$p.value
)
names(pValues) <- names(medianRatios)

names(pValues) <- names(medianRatios)
qValues <- data.frame(apply(pValues,2,p.adjust,method="BH"))	

sqa <- list()
class(sqa) <- "safeQuantAnalysis"

sqa$eset <- esetProteinRatios
sqa$qValue <- qValues
sqa$pValue <- pValues
sqa$ratio <- medianRatios

### GRAPIHCS

pdf(pdfFile)

par(mfrow = c(2,2))
plotExpDesign(esetR1)
plotExpDesign(esetR2)
plotExpDesign(esetR3)
plotExpDesign(esetR4)

par(mfrow=c(1,1))

venn(list(r1 = as.character(fData(sqaProtR1$eset)$proteinName)
				,r2 = as.character(fData(sqaProtR2$eset)$proteinName)
				,r3 = as.character(fData(sqaProtR3$eset)$proteinName)
				,r4 = as.character(fData(sqaProtR4$eset)$proteinName)
				)
)
mtext(side=3,text="Protein")

venn(list(r1 = unique(as.character(fData(esetR1)$peptide))
				,r2 = unique(as.character(fData(esetR2)$peptide))
				,r3 = unique(as.character(fData(esetR3)$peptide))
				,r4 = unique(as.character(fData(esetR4)$peptide))
		
		)
)
mtext(side=3,text="Peptide")

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

par(mfrow=c(2,2))
.allpValueHist(sqa)

par(mfrow=c(2,2))
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=T, isLegend=T ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=T, isLegend=F, pvalCutOff=qValueThrs)

plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=F , isLegend=F ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=F, isLegend=F ,pvalCutOff=qValueThrs)

# plot volcanoes
if(T){
	# volcanoes qvalue
	par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.3 )
	for(i in 1:ncol(qValues)){
		caseCondition <- names(medianRatios)[i]		
		plot(medianRatios[,i],-log10(qValues[,i])
				, pch=19
				, xlab= paste("log2(",caseCondition,"/WT",")", sep="" )
				, ylab= paste("-log10(qValue)",sep="")
				,col=c("black","blue")[(qValues[,i] < qValueThrs)+1 ]
		)
		grid()
	}
}

dev.off()

### GRAPIHCS END

### XLS EXPORT

ratios <- exprs(esetProteinRatios)
colnames(ratios) <- paste("log2ratio",colnames(exprs(esetProteinRatios)),sep="_")
names(medianRatios) <-  paste("median_log2ratio",names(medianRatios),sep="_")
names(pValues) <-  paste("pValue",names(pValues),sep="_")
names(qValues) <-  paste("qValue",names(pValues),sep="_")
out <- cbind(fData(esetProteinRatios),ratios,medianRatios,pValues,qValues)
write.table(file=xlsFile,out,row.names=F,sep="\t")

### XLS EXPORT END



