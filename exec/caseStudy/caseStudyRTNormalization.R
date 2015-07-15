# Evaluate Accuracy and Precision when applying "Global Median" and Retention Time Dependent Normalization.
# 
# Author: erikahrne
###############################################################################

plotBarheVolcano <- function(sqa, cond=1){
	
	plot(sqa$ratio[,cond],-log10(sqa$qValue[,cond]), pch=19, col=fData(sqa$eset)$isNormAnchor+1, ylab="-log10(qValue)", xlab="log2(ratio)", xlim=c(0,max(sqa$ratio[,cond]*1.1)),xaxs="i")
	legend("bottomright",c("HUMAN","BARHE"),pch=19,col=c(2,1), cex=1.5)
	
}


plotBarheRoc <- function(sqa, pValueCutOffs=c(seq(0,0.1,length=1000),seq(0.1,1,length=100)), cond=1,type="l"){
	
	isBarhe <- !fData(sqa$eset)$isNormAnchor
	
	tpr <- c() # TP/(TP+FN) # a.k.a sensitivity
	fpr <- c() # FP/(FP+TN) # a.k.a 1 - specificity
	
	for(pc in pValueCutOffs){
		
		tp <- sum((sqa$pValue[,cond] < pc) & isBarhe)
		fp <- sum((sqa$pValue[,cond] < pc) & !isBarhe)
		tn <- sum((sqa$pValue[,cond] > pc) & !isBarhe)
		fn <- sum((sqa$pValue[,cond] > pc) & isBarhe)
		
		#tpr <- c(tpr,tp/(tp+fn))
		tpr <- c(tpr,tp/sum(isBarhe))
		#fpr <- c(fpr,fp/(fp+tn))
		fpr <- c(fpr,fp/sum(!isBarhe))
		
	}
	
	plot(fpr,tpr,type=type, lwd=1.5, main=names(sqa$pValue)[cond])
	
	out <- list()
	
	out$pValueCutOffs <- pValueCutOffs 
	out$fpr <- fpr 
	out$tpr <- tpr 
	
	return(out)
	
}


source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/caseStudy/initCaseStudySession.R")


### PROTEOME MIX TEST CASE
### @TODO evaluate new normalization strategy  
# In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. 2013
# normalize for variations in elelctrospray ionization current
# plot ratio as a function of retention time

### INIT



#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/old/Peptides_All_Max_Sens_241011_FILTERED.csv"
#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/old/peptide_FILTERED.csv"
#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/old/peptides1_FILTERED.csv"
#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/old/peptides2_FILTERED.csv"

#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/HEK-Nucleases-110814-Peptides_FILTERED.csv"

#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/Cecile-CS-Yewtivya-MSA-Features_FILTERED.csv"
#progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/Cecile-CS-Yewtivya-HCD-Features_FILTERED.csv"

progenesisPeptideCsvFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/peptides_proteomeMix_FILTERED.csv"

fdrCutOff <- 0.01
qValueThrs <-  0.05
pMassTolWindow <-  c(0,8)

#pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/tests/caseStudy/pdf/rtNormalization.pdf"
pdfFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/tests/caseStudy/pdf/rtNormalization_vsn.pdf"

### INIT END

eset <- parseProgenesisPeptideCsv(file=progenesisPeptideCsvFile,expDesign= getExpDesignProgenesisCsv(progenesisPeptideCsvFile))
fData(eset)$isNormAnchor <- grepl("HUMAN",fData(eset)$proteinName) 

#filter
eset <- eset[(fData(addIdQvalues(eset))$idQValue < fdrCutOff)
				& !isDecoy(fData(eset)$proteinName) 
				& !isCon(fData(eset)$proteinName)
				& ((fData(eset)$pMassError > pMassTolWindow[1]) &&  (fData(eset)$pMassError < pMassTolWindow[2]))
				, ] # id score

		
plotPrecMassErrorDistrib(fData(eset)$pMassError, isDecoy(fData(eset)$proteinName),pMassTolWindow=pMassTolWindow)		
#plotPrecMassErrorVsScore(fData(eset)$pMassError, fData(eset)$idScore, isDecoy(fData(eset)$proteinName), pMassTolWindow=pMassTolWindow )		
		
### rt normalize
rtNormFactors <- getRTNormFactors(eset)
esetNorm <- rtNormalize(eset,rtNormFactors)

### @TODO TMP
#library(vsn)
#esetNorm <- justvsn(eset)
#exprs(esetNorm) <- 2^exprs(esetNorm)

# global median normalization
eset <- normalize(eset, method="global")



### roll up
#esetPeptide <- rollUp(eset[!fData(eset)$isNormAnchor,],featureDataColumnName= c("proteinName"),isProgressBar=T )
#esetPeptideNorm <- rollUp(esetNorm[!fData(esetNorm)$isNormAnchor,],featureDataColumnName= c("proteinName"),isProgressBar=T)

esetPeptide <- rollUp(eset,featureDataColumnName= c("proteinName"),isProgressBar=T )
esetPeptideNorm <- rollUp(esetNorm,featureDataColumnName= c("proteinName"),isProgressBar=T)


### stat analysis
sqaPeptide	<- safeQuantAnalysis(esetPeptide, method=c("naRep"))
sqaPeptideNorm	<- safeQuantAnalysis(esetPeptideNorm, method=c("naRep") )

### GRAPHICS

pdf(pdfFile)

hist(round(fData(eset)$retentionTime), breaks=length(unique(round(fData(eset)$retentionTime))))
abline(h=100)

par(mfrow=c(1,1))

plotRTNormSummary(rtNormFactors, lwd=1.5, col=as.character(.getConditionColors(eset)[pData(eset)[names(rtNormFactors),]$condition,]),condNames=as.character(unique(pData(eset)$condition)) )
plotRTNorm(rtNormFactors,eset, ylim=c(-4,4),samples=2)

par(mfrow=c(1,2), mar=c(6.1,4.1,4.1,2.1), las=2)
.qcPlots(eset[fData(eset)$isNormAnchor,], selection=2)
.qcPlots(esetNorm[fData(esetNorm)$isNormAnchor,], selection=2)

boxplot(sqaPeptide$cv, ylim=c(0,0.5))
boxplot(sqaPeptideNorm$cv, ylim=c(0,0.5))

boxplot(sqaPeptide$cv- sqaPeptideNorm$cv, ylim=c(-0.3,0.3), ylab="Delta C.V.", las=2)
boxplot(sqaPeptide$cv[!fData(sqaPeptide$eset)$isNormAnchor,]- sqaPeptideNorm$cv[!fData(sqaPeptideNorm$eset)$isNormAnchor,], ylim=c(-0.3,0.3), ylab="Delta C.V.", las=2)
abline(h=0)

### Valid D.E. Peptides


nbDePeptides <- rbind(apply(sqaPeptide$qValue[!fData(sqaPeptide$eset)$isNormAnchor,] < qValueThrs,2,sum,na.rm=T),apply(sqaPeptideNorm$qValue[!fData(sqaPeptideNorm$eset)$isNormAnchor,] < qValueThrs,2,sum,na.rm=T))

barplot(nbDePeptides,beside=T, col=3:4, las=2)
legend("left",c("old Norm","new norm"), fill=3:4)

par(mfrow=c(1,2), mar=c(12.1,4.1,4.1,2.1), las=2)
boxplot(sqaPeptide$ratio[!fData(sqaPeptide$eset)$isNormAnchor,], ylim=c(-5,1))
abline(h=c(log2(0.02/0.1) , log2(0.015/0.1), log2(0.01/0.1)))
boxplot(sqaPeptideNorm$ratio[!fData(sqaPeptideNorm$eset)$isNormAnchor,], ylim=c(-5,1))
abline(h=c(log2(0.02/0.1) , log2(0.015/0.1), log2(0.01/0.1)))


### ROC curve
for(i in 1:ncol(sqaPeptide$ratio)){
	
	out <- plotBarheRoc(sqaPeptide, cond=i, ,type="n")
	outNorm <- plotBarheRoc(sqaPeptideNorm, cond=i, ,type="n")
	
	plot(out$fpr,out$tpr,type="l", col="red", main=names(sqaPeptide$ratio)[i])
	lines(outNorm$fpr,outNorm$tpr, col="blue")
	legend("bottomright",c("global norm","rt norm"),lty=1, col=c("red","blue"))
}

print("DONE")

dev.off()

### GRAPHICS OFF
if(F){
	library(vsn)
	
	
	esetVsnNorm <- justvsn(eset)
	plot(log2(exprs(eset[,1])), exprs(esetVsnNorm[,1]))
	
	plot(exprs(eset[,1]), exprs(eset[,2]),log="xy")
	plot(exprs(esetVsnNorm[,1]), exprs(esetVsnNorm[,2]) ,log="xy")
	
	data("kidney")
	
	
	plot(exprs(kidney)[,1],exprs(kidney)[,2], log="xy")
	
	kidneyVSN <- justvsn(kidney)
	plot(exprs(kidneyVSN)[,1],exprs(kidneyVSN)[,2], log="xy")
	
	plot(log2(exprs(kidney[,1])), exprs(kidneyVSN[,1]))
	
}
