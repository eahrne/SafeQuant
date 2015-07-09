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

### PARAMS

progenesisFileRep1 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Phospho/Phospho_Set_1.csv"
progenesisFileRep2 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Phospho/Phospho_Set_2.csv"
progenesisFileRep3 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Phospho/Phospho_Set_3_noGroupedAccessions_4.csv"
progenesisFileRep4 <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/MRuegg/ChristianZimmerli_111/AdditionalFiles/Pankaj_Phospho/Phospho_Set_4.csv"

fastaFile <- "/Volumes/pcf01$/Schmidt_Group/Databases/Latest_UniProt/uniprot-Mouse_301014.fasta"

rDataFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/data/caseStudyZimmerliMergeAndQunatifyPO4.rData"

#pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyPO4/sq_intersect.pdf"
#xlsFile <-  "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyPO4/sq_intersect.xls"

pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyPO4/sq.pdf"
xlsFile <-  "/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/out/caseStudyZimmerliMergeAndQuantifyPO4/sq.xls"

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
	
	esetPeptideR1 <- rollUp(esetR1,featureDataColumnName =  c("peptide","ptm"))
	esetPeptideR2 <- rollUp(esetR2,featureDataColumnName =  c("peptide","ptm"))
	esetPeptideR3 <- rollUp(esetR3,featureDataColumnName =  c("peptide","ptm"))
	esetPeptideR4 <- rollUp(esetR4,featureDataColumnName =  c("peptide","ptm"))
	
	### ROLL-UP END
	
	### SQA
	
	sqaPeptideR1 <- safeQuantAnalysis(esetPeptideR1 )
	sqaPeptideR2 <- safeQuantAnalysis(esetPeptideR2 )
	sqaPeptideR3 <- safeQuantAnalysis(esetPeptideR3)
	sqaPeptideR4 <- safeQuantAnalysis(esetPeptideR4)
	
	### SQA END 
		
	proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
	# dirty fix check if ACs in Progenesis file are stripped
	names(proteinDB) <- stripACs(names(proteinDB))
 	
	save(esetR1,esetR2,esetR3,esetR4,esetPeptideR1,esetPeptideR2,esetPeptideR3,esetPeptideR4,sqaPeptideR1,sqaPeptideR2,sqaPeptideR3,sqaPeptideR4,proteinDB,file=rDataFile)
	#save(esetR1,esetR2,esetR4,esetPeptideR1,esetPeptideR2,esetPeptideR4,sqaPeptideR1,sqaPeptideR2,sqaPeptideR4,file=rDataFile)
	
}else{
	
	load(rDataFile)
	
}


allPeptides <- unique(c( rownames(sqaPeptideR1$eset)
						,rownames(sqaPeptideR2$eset)
						,rownames(sqaPeptideR3$eset)
						,rownames(sqaPeptideR4$eset)
					))

		
intersectPeptides <- intersect(
						intersect(rownames(sqaPeptideR1$eset),rownames(sqaPeptideR2$eset))
						,intersect(rownames(sqaPeptideR3$eset),rownames(sqaPeptideR4$eset))
				)

#intersectPeptides <- intersect(
#				intersect(rownames(sqaPeptideR1$eset),rownames(sqaPeptideR2$eset))
#				,rownames(sqaPeptideR4$eset))
			
				
#selPeptides	<- 	intersectPeptides		
selPeptides	<- 	allPeptides		

expMatrixRatios <- data.frame(sqaPeptideR1$ratio[selPeptides,]
			,sqaPeptideR2$ratio[selPeptides,]
			,sqaPeptideR3$ratio[selPeptides,]
			,sqaPeptideR4$ratio[selPeptides,])
rownames(expMatrixRatios) <- selPeptides

expMatrixRatios <- expMatrixRatios[,rev(order(names(expMatrixRatios)))]

#head(expMatrixRatios)
#head(expMatrixRatios[,rev(order(names(expMatrixRatios)))])

featureData <- unique(rbind(fData(sqaPeptideR1$eset)[,c(1:3,9)]
				,fData(sqaPeptideR2$eset)[,c(1:3,9)]
				,fData(sqaPeptideR3$eset)[,c(1:3,9)]
				,fData(sqaPeptideR4$eset)[,c(1:3,9)]
		
		))[selPeptides,]
rownames(featureData) <- selPeptides



expDesign <- data.frame(row.names=names(expMatrixRatios),isControl=rep(F,ncol(expMatrixRatios)),condition=gsub("\\_[0-9]{1,}","",names(expMatrixRatios)) )
esetPeptideRatios <- createExpressionDataset(expressionMatrix=as.matrix(expMatrixRatios),expDesign=expDesign,featureAnnotations=featureData)		 

### add motif-x 
esetPeptideRatios <- .addPTMCoord(esetPeptideRatios,proteinDB,motifLength=4, isProgressBar=T)

#pData(esetPeptideRatios)


# caluc median ratios per condition (merged data)
medianRatios <- getSignalPerCondition(esetPeptideRatios)

# p-values
pValues <- data.frame(eBayes(lmFit(esetPeptideRatios[,1:4 ]))$p.value
		,eBayes(lmFit(esetPeptideRatios[,5:8 ]))$p.value
		,eBayes(lmFit(esetPeptideRatios[,9:12 ]))$p.value
)
#pValues <- data.frame(eBayes(lmFit(esetPeptideRatios[,1:3 ]))$p.value
#		,eBayes(lmFit(esetPeptideRatios[,4:6 ]))$p.value
#		,eBayes(lmFit(esetPeptideRatios[,7:9 ]))$p.value
#)

names(pValues) <- names(medianRatios)

names(pValues) <- names(medianRatios)
qValues <- data.frame(apply(pValues,2,p.adjust,method="BH"))	

sqa <- list()
class(sqa) <- "safeQuantAnalysis"

sqa$eset <- esetPeptideRatios
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


venn(list(r1 = unique(rownames(esetPeptideR1))
				,r2 = unique(rownames(esetPeptideR2))
				,r3 = unique(rownames(esetPeptideR3))
				,r4 = unique(rownames(esetPeptideR4))
		
		)
)
mtext(side=3,text="PO4 Peptide")

conditionColors =.getConditionColors(esetPeptideRatios)
samplecolors =  as.vector(unlist(conditionColors[pData(esetPeptideRatios)$condition,]))
breaks=seq(-1.5,1.5,length=20)

m <- exprs(esetPeptideRatios)
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
legend("left",levels(pData(esetPeptideRatios)$condition), fill=as.character(conditionColors[,1]), cex=0.7, box.col=0)

par(mfrow=c(2,2))
.allpValueHist(sqa)

par(mfrow=c(2,2))
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=T, isLegend=T ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=T, isLegend=F, pvalCutOff=qValueThrs)

plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="UP", isAdjusted=F , isLegend=F ,pvalCutOff=qValueThrs)
plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,1), main="DOWN", isAdjusted=F, isLegend=F ,pvalCutOff=qValueThrs)

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

ratios <- exprs(esetPeptideRatios)
colnames(ratios) <- paste("log2ratio",colnames(exprs(esetPeptideRatios)),sep="_")
names(medianRatios) <-  paste("median_log2ratio",names(medianRatios),sep="_")
names(pValues) <-  paste("pValue",names(pValues),sep="_")
names(qValues) <-  paste("qValue",names(pValues),sep="_")
out <- cbind(fData(esetPeptideRatios),ratios,medianRatios,pValues,qValues)
write.table(file=xlsFile,out,row.names=F,sep="\t")

### XLS EXPORT END




