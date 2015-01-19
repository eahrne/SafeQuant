#!/usr/bin/Rscript

#  Run Safequant
# Identification FDR calulcation: Calculate fdr and filter before moving on to the next level, spectrum, peptide, protein 
# If not large proteins are very much favored
# TEST
#"/Volumes/pcf01\$/Schmidt_Group/ProjectSQ/DBumann/DanielaMariaRemus_5/20130711-143152_DR1/Progenesis/Bumann_DMR_rpoE-phoP-mCherry_191113/proteins.csv"
#"/Volumes/pcf01\$/Schmidt_Group/ProjectSQ/DBumann/DanielaMariaRemus_5/20130711-143152_DR1/Progenesis/Bumann_DMR_rpoE-phoP-mCherry_191113/peptides.csv"

# Author: erikahrne
###############################################################################


# TEST FILE
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/peptides1_FILTERED.csv /Volumes/pcf01\$/Schmidt_Group/Databases/SwissProt_Databases/s_human_d_201405.fasta
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/proteins1.csv
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/TMT_6-Plex_Scaffold_Raw_Export_Example.xls

############################################################### INIT ############################################################### 

if(F){
	### 
	library("SafeQuant")
}else if(file.exists("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")){
	#### DEPENDANCIES
	
	#@TEMP
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")
}else{
	#@TEMP TPP
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/UserOptions.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/ExpressionAnalysis.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/SafeQuantAnalysis.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/Graphics.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/IdentificationAnalysis.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/Parser.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/TMT.R")
	source("/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/UserOptions.R")
	
}

VERSION <- 2.01

### USER CMD LINE OPTIONS
userOptions <- getUserOptions(version=VERSION)
### USER CMD LINE OPTIONS END

if(userOptions$verbose) print(userOptions$proteinQuant)

### SUPRESS WARNINGS
if(!userOptions$verbose){
	options(warn=-1)
}

#### DEPENDENCIES

if(userOptions$verbose) print(userOptions)

if(userOptions$verbose) cat("LOADING DEPENDENCIES \n")	

suppressPackageStartupMessages(library("affy", quiet=T))
suppressPackageStartupMessages(library("limma", quiet=T))
suppressPackageStartupMessages(library(gplots, quiet=T)) # volcano plot
suppressPackageStartupMessages(library(seqinr, quiet=T))
suppressPackageStartupMessages(library(corrplot, quiet=T))

############################################################### PARSING ############################################################### 

if(userOptions$verbose) cat("PARSING INPUT FILE \n")	

# get file type
fileType <- .getFileType(userOptions$inputFile)

normMethod <- c("global")

### Progenesis Export
if(fileType %in% c("ProgenesisProtein","ProgenesisFeature")){
	
	# default
	expDesign <- getExpDesignProgenesisCsv(userOptions$inputFile)
	
	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)		
	}

	if(fileType == "ProgenesisProtein"){
		cat("INFO: PARSING PROGENESIS PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisProteinCsv(file=userOptions$inputFile,expDesign=expDesign)
		
	}else{ #"ProgenesisFeature"
		cat("INFO: PARSING PROGENESIS FEATURE EXPORT FILE ",userOptions$inputFile, "\n" )
		normMethod <- c("rt")
		eset <- parseProgenesisFeatureCsv(file=userOptions$inputFile,expDesign=expDesign)
	}
	
# Scaffold Export (TMT data)
}else if(fileType == "ScaffoldTMT"){
	cat("INFO: PARSING SCAFFOLD RAW EXPORT FILE ",userOptions$inputFile, "\n" )
	
	# get default experimental design
	# six plex or ten plex ?
	# use default experimental design unless specified by the user
	if(.getNbPlex(userOptions$inputFile) == 6){
		# 6-plex default: 1,2,3:4,5,6 
		expDesign <- data.frame(condition=paste("Condition",sort(rep(c(1,2),3))),isControl=sort(rep(c(T,F),3),decreasing=T) )
	}else{
		# 10-plex default is "1,4,7,10:2,5,8:3,6,9"
		expDesign <- data.frame(condition=paste("Condition",c(1,2,3,1,2,3,1,2,3,1)),isControl=c(T,F,F,T,F,F,T,F,F,T) )
	}
	
	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)
	}
	
	eset <- parseScaffoldRawFile(file=userOptions$inputFile,expDesign=expDesign)
}else if(fileType == "GenericCSV"){
	
}else{
	stop("Unknown File Type", userOptions$inputFile)
}

# test option, limit number of entries
if(userOptions$test){
	eset <- eset[1:nrow(eset) %in% sample(1:min(300,c(nrow(eset))),replace=F),]
}

# parse .fasta file
if(!is.na(userOptions$proteinFastaFile)){
	cat("INFO: PARSING PROTEIN SEQUENCE DB ",userOptions$proteinFastaFile, "\n" )
	### read protein db
	proteinDB <- read.fasta(userOptions$proteinFastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
}

############################################################### CREATE DATA MODEL ############################################################### 

if(userOptions$verbose) print(eset)
if(userOptions$verbose) print(pData(eset))
if(userOptions$verbose) print(names(fData(eset)))

#### CREATE FEATURE DATA AND FILTER (pre-rollup)

# generic
filter <- data.frame(
		con=isCon(fData(eset)$proteinName)	# contaminants	
		,ac = !(grepl(userOptions$selectedProteinName,fData(eset)$proteinName,ignore.case=T)) # protein ac
)

if("pMassError" %in% names(fData(eset))){
	### applicable to Progenesis feature Exports	
	
	filter <- cbind(filter, pMassError=
					(fData(eset)$pMassError < userOptions$precursorMassFilter[1])
					| (fData(eset)$pMassError > userOptions$precursorMassFilter[2]) # precursor mass tolerance
	)
}

if("ptm" %in% names(fData(eset))){
	
	# add motif-X and ptm coordinates
	if(exists("proteinDB")){
		cat("INFO: EXTRACTING PTM COORDINATES AND MOTIFS\n")
		eset <- .addPTMCoord(eset,proteinDB,motifLength=4, isProgressBar=T)
	}
	filter <- cbind(filter
					, ptm = !(grepl(userOptions$selectedModifName,as.character(fData(eset)$ptm),ignore.case=T))
					, nbPtmsPerPeptide = (fData(eset)$nbPtmsPerPeptide > userOptions$maxNbPtmsPerPeptide) )

}

if(!("nbPeptides" %in% names(fData(eset)))){
	### set nb peptides per protein
	eset <- setNbPeptidesPerProtein(eset)
}
filter <- cbind(filter,nbPeptides=(fData(eset)$nbPeptides < userOptions$minNbPeptidesPerProt))

if("idScore" %in% names(fData(eset))){
	eset <- addIdQvalues(eset)
	filter <- cbind(filter,qvalue=fData(eset)$idQValue > userOptions$fdrCutoff)
}

# set pre-rollup filters
eset <- .setFilter(eset,filter=filter)

### make sure at least 1 feature pass the filter
if(sum(!fData(eset)$isFiltered,na.rm=T) == 0){
	stop("CHECK FILTER SETTINGS. ALL FEATURES WERE FILTERED OUT")
}

#### CREATE FEATURE DATA AND FILTER END

### SET ANCHOR PROTEINS
fData(eset)$isNormAnchor <- grepl(userOptions$normAC,fData(eset)$proteinName)

if(userOptions$verbose){
	cat("\nNB. ANCHOR PROTEINS: ")
	cat(sum(fData(eset)$isNormAnchor))
	cat("\n")
	print(fData(eset)$proteinName[fData(eset)$isNormAnchor])
	cat("\n")
}

### SET ANCHOR PROTEINS END

############################################################### EXPRESSION ANALYSIS ############################################################### 

### non-pairwise stat test
statMethod <- c("")
if(userOptions$SNonPairWiseStatTest) statMethod <- c("all")

### normalize eset
#cat("INFO: NORMALIZING DATA USING '",normMethod ,"' METHOD \n",sep="")

### if user has specified anchor proteins -> always use global method 
if(sum(fData(eset)$isNormAnchor == FALSE, na.rm=T ) > 0 ){
	cat("INFO: 'GLOBAL' Normalization always used in combination with SAnchorProtein option \n")	
	normMethod <- "global"
}
esetNorm <- normalize(eset,method=normMethod)
### add pseudo (baseline) intensity
baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(esetNorm)[,1])),promille=5)
exprs(esetNorm)[is.na(exprs(esetNorm)) | (exprs(esetNorm) < 0)  ] <- 0 
exprs(esetNorm) <- exprs(esetNorm) + baselineIntensity

#exprs(esetNorm)[is.na(exprs(esetNorm)) | (exprs(esetNorm) < 0)  ] <- 0 
#baselineIntensity <- .getBaselineIntensityTMT(exprs(esetNorm)[,1]) 
#exprs(esetNorm)[is.na(exprs(esetNorm)) | (exprs(esetNorm) < baselineIntensity)  ] <- baselineIntensity 

### ProgenesisProtein -> NO ROLL-UP
### ProgenesisFeature -> ROLL-UP PEPTIDE (ROLL-UP PROTEIN UNLESS USER SPECIFIED )
### ScaffoldTMT -> ROLL-UP PROTEIN

if((fileType == "ProgenesisProtein")){
	
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered  | isDecoy(fData(esetNorm)$proteinName)
	sqaProtein <- safeQuantAnalysis(esetNorm, method=statMethod)
}else if((fileType == "ScaffoldTMT")){

	# roll-up protein level
	cat("INFO: ROLL-UP PROTEIN LEVEL\n")
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered |  isDecoy(fData(esetNorm)$proteinName)
	sqaProtein <- safeQuantAnalysis(rollUp(esetNorm,featureDataColumnName= c("proteinName"),isProgressBar=T), method=statMethod)
	fData(sqaProtein$eset)$isFiltered <- fData(sqaProtein$eset)$isFiltered | isDecoy(fData(sqaProtein$eset)$proteinName)
	
}else{
	
	# roll-up peptide level
	cat("INFO: ROLL-UP PEPTIDE LEVEL\n")
	esetPeptide <- rollUp(esetNorm,featureDataColumnName= c("peptide","ptm"),isProgressBar=T)
	# fdr filter
	# replace qValues by rollUp level qValues ()
	esetPeptide <- addIdQvalues(esetPeptide)
	# update filter to exclude peptide level hight qValues
	fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | (fData(esetPeptide)$idQValue > userOptions$fdrCutoff)
	
	if(userOptions$proteinQuant){
		cat("INFO: ROLL-UP PROTEIN LEVEL\n")
		esetProtein <- rollUp(esetPeptide,featureDataColumnName= c("proteinName"),isProgressBar=T)
		esetProtein <- addIdQvalues(esetProtein)
		fData(esetProtein)$isFiltered <- fData(esetProtein)$isFiltered | (fData(esetProtein)$idQValue > userOptions$fdrCutoff) | isDecoy(fData(esetProtein)$proteinName)
		sqaProtein <- safeQuantAnalysis(esetProtein, method=statMethod)

	}
	
	fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | isDecoy(fData(esetPeptide)$proteinName)
	sqaPeptide <- safeQuantAnalysis(esetPeptide, method=statMethod)
	
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered | isDecoy(fData(esetNorm)$proteinName)
	
	if(userOptions$top3){
		cat("INFO: ROLL-UP TOP3\n")
		esetTop3 <-  rollUp(esetPeptide,featureDataColumnName= c("proteinName"),isProgressBar=T, method="top3")
	}
}

### IBAQ
if(userOptions$iBAQ){
	cat("INFO: CALCULATING IBAQ VALUES\n")
	if(exists("proteinDB")){
		esetIBAQ <-  getIBAQEset(sqaProtein$eset, proteinDB=proteinDB)
	}else{
		cat("ERROR: proteinDB NOT FOUND NO iBAQ VALUES CALCULATED\n")
	}
}

### EXPRESSION ANALYSIS END

############################################################### EXPORTS ############################################################### 

#### SET WORKING DIR

if(!file.exists(userOptions$outputDir)) dir.create(userOptions$outputDir)
if(userOptions$verbose) cat("INFO: CREATED DIRECTORY",  userOptions$outputDir,"\n")
#### SET WORKING DIR

##I/O: set export file paths
userOptions$pdfFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,".pdf",sep=""))
userOptions$peptideTsvFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_PEPTIDE.tsv",sep=""))
userOptions$proteinTsvFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_PROTEIN.tsv",sep=""))
userOptions$paramsFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_SQ_PARAMS.TXT",sep=""))
userOptions$rDataFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_SQ.rData",sep=""))
### GRAPHICS

pdf(userOptions$pdfFile)
parDefault <- par()
CONDITIONCOLORS <- .getConditionColors(esetNorm)

### EXPDESIGN PLOT
plotExpDesign(esetNorm, version=VERSION)
### EXPDESIGN PLOT END

# l <- layout(rbind(c(1,1), c(2,3)))
# layout.show(l)

### IDENTIFICATION PLOTS
#if(fileType == "ProgenesisProtein") layout(rbind(c(1, 1), c(2, 3)) )
if(fileType == "ProgenesisProtein") par(mfrow=c(2,2))
if(fileType == "ScaffoldTMT") par(mfrow=c(2,2))
#if(fileType == "ProgenesisFeature")layout(rbind(c(1,1,1,2,2,2), c(3,3, 4,4,5,5)))
if(fileType == "ProgenesisFeature") par(mfrow=c(2,3))
.idOverviewPlots()
if(fileType == "ProgenesisFeature")par(mfrow=c(2,2))
if(exists("sqaPeptide")) .idPlots(sqaPeptide$eset, selection=c(1,3), main="Peptide Level", qvalueThrs=userOptions$fdrCutoff)
if(exists("sqaProtein")) .idPlots(sqaProtein$eset, selection=c(1,3), main="Protein Level", qvalueThrs=userOptions$fdrCutoff)
par(parDefault)
### IDENTIFICATIONS PLOTS END

### QUANT. QC PLOTS 

rowSelEsetNorm <- sample(nrow(esetNorm),min(c(500,nrow(esetNorm))) ,replace=F)
rowSelEset <- sample(nrow(eset),min(c(1000,nrow(eset))) ,replace=F)

### CORRELATION PLOTS
### COR OR PAIRS PLOT. IF FEWER THAN X SAMPLES 
if(ncol(esetNorm) < 8){
	pairsAnnot(log10(exprs(esetNorm))[rowSelEsetNorm,],textCol=as.character(CONDITIONCOLORS[pData(esetNorm)$condition,]))
}else{
	.correlationPlot(log10(exprs(esetNorm))[rowSelEsetNorm,], labels=as.character(unique(pData(esetNorm)$condition)), textCol=as.character(CONDITIONCOLORS[pData(esetNorm)$condition,]))
}
### COR OR PAIRS PLOT. IF FEWER THAN X CONDITIONS 
if(length(unique(pData(esetNorm)$condition)) < 8){
	pairsAnnot(log10(getSignalPerCondition(esetNorm[rowSelEsetNorm,]))[,unique(pData(esetNorm)$condition) ],textCol=as.character(CONDITIONCOLORS[unique(pData(esetNorm)$condition),]))
}else{
	.correlationPlot(log10(getSignalPerCondition(esetNorm[rowSelEsetNorm,]))[,unique(pData(esetNorm)$condition) ],textCol=as.character(CONDITIONCOLORS[unique(pData(esetNorm)$condition),]))
}

par(parDefault)

if(fileType == "ProgenesisFeature"){
	
	par(mfrow=c(2,1), mar=c(4.5,6.1,4.1,6.1))
	plotPrecMassErrorDistrib(eset, pMassTolWindow=userOptions$precursorMassFilter)
	plotPrecMassErrorVsScore(eset[rowSelEset,], pMassTolWindow=userOptions$precursorMassFilter)
	par(parDefault)
}

layout(rbind(c(1,2), c(3,3)))
missinValueBarplot(eset)
barplotMSSignal(eset)
par( mar=c(6.5,5.1,2.5,3.1))
cvBoxplot(esetNorm)
par(parDefault)

### retention time plot
if("rt" %in% normMethod ) plotRTNormSummary(eset,lwd=2)

### QUANT. QC PLOTS END

par(mfrow=c(1,2))
#if(userOptions$eBayes | (exists("sqaPeptide") & exists("sqaProtein"))) par(mfrow=c(2,2))

### QUANT. STAT. PLOTS 

### VAILD FEATURES VS. pValue/qValue

if(exists("sqaProtein")){
	
	plotNbValidDeFeaturesPerFDR(sqaProtein,
			upRegulated=F
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=T
			,isAdjusted=T
			,ylab="Protein Counts"
			,main="DOWN REGULATION"
	)
	
	plotNbValidDeFeaturesPerFDR(sqaProtein,
			upRegulated=T
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=F
			,isAdjusted=T
			,ylab="Protein Counts"
			,main="UP REGULATION"
	)

}else if(exists("sqaPeptide")){
	
	plotNbValidDeFeaturesPerFDR(sqaPeptide,
			upRegulated=F
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=T
			,isAdjusted=T
			,ylab="Peptide Counts"
			,main="DOWN REGULATION"
	)
	
	plotNbValidDeFeaturesPerFDR(sqaPeptide,
			upRegulated=T
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=F
			,isAdjusted=T
			,ylab="Peptide Counts"
			,main="UP REGULATION"
	)
	
}

if(userOptions$eBayes){ 
	if(exists("sqaProtein")){
		
		plotNbValidDeFeaturesPerFDR(sqaProtein,
				upRegulated=F
				,log2RatioCufOff=log2(userOptions$ratioCutOff)
				,pvalRange=c(0,0.15)
				,pvalCutOff=userOptions$deFdrCutoff
				,isLegend=T
				,isAdjusted=F
				,ylab="Protein Counts"
				,main="DOWN REGULATION"
		)
		
		plotNbValidDeFeaturesPerFDR(sqaProtein,
				upRegulated=T
				,log2RatioCufOff=log2(userOptions$ratioCutOff)
				,pvalRange=c(0,0.15)
				,pvalCutOff=userOptions$deFdrCutoff
				,isLegend=F
				,isAdjusted=F
				,ylab="Protein Counts"
				,main="UP REGULATION"
		)
	}else if(exists("sqaPeptide")){
		
		plotNbValidDeFeaturesPerFDR(sqaPeptide,
				upRegulated=F
				,log2RatioCufOff=log2(userOptions$ratioCutOff)
				,pvalRange=c(0,0.15)
				,pvalCutOff=userOptions$deFdrCutoff
				,isLegend=T
				,isAdjusted=F
				,ylab="Peptide Counts"
				,main="DOWN REGULATION"
		)
		
		plotNbValidDeFeaturesPerFDR(sqaPeptide,
				upRegulated=T
				,log2RatioCufOff=log2(userOptions$ratioCutOff)
				,pvalRange=c(0,0.15)
				,pvalCutOff=userOptions$deFdrCutoff
				,isLegend=F
				,isAdjusted=F
				,ylab="Peptide Counts"
				,main="UP REGULATION"
		)
	}
} 

par(parDefault)

if(exists("sqaProtein")){
	hClustHeatMap(sqaProtein$eset,main="Protein Level")
	
	plotVolcano(sqaProtein
			, main="Protein Level" 
			, ratioThrs= userOptions$ratioCutOff
			, pValueThreshold= userOptions$deFdrCutoff
			, adjusted = T)
}else if(exists("sqaPeptide")){
	hClustHeatMap(sqaPeptide$eset,main="Peptide Level")
	
	plotVolcano(sqaPeptide
			, main="Peptide Level" 
			, ratioThrs= userOptions$ratioCutOff
			, pValueThreshold= userOptions$deFdrCutoff
			, adjusted = T)
} 

if(userOptions$eBayes){
	
	if(exists("sqaProtein")){
		
		plotVolcano(sqaProtein
				, main="Protein Level"
				, ratioThrs= userOptions$ratioCutOff
				, pValueThreshold= userOptions$deFdrCutoff
				, adjusted = F)
		
		par(mfrow=c(2,2))
		if(nrow(CONDITIONCOLORS) > 4) par(mfrow=c(3,3))
		.allpValueHist(sqaProtein)
		par(parDefault)
		
	}else  if(exists("sqaPeptide")){
		par(mfrow=c(2,2))
		if(nrow(CONDITIONCOLORS) > 4) par(mfrow=c(3,3))
		.allpValueHist(sqaPeptide)
		par(parDefault)
		
		plotVolcano(sqaPeptide
					, main="Peptide Level"
					, ratioThrs= userOptions$ratioCutOff
					, pValueThreshold= userOptions$deFdrCutoff
					, adjusted = F)
	}
}

### QUANT. STAT. PLOTS END

par(mfrow=c(1,1))

cat("INFO: CREATED FILE ", userOptions$pdfFile,"\n")

graphics.off()

### GRAPHICS END

### TSV EXPORT

if(exists("sqaPeptide")){ 
	
	selFDataCol <- c("peptide","proteinName","proteinDescription", "idScore","idQValue"
						,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures", "isNormAnchor") 
	selFDataCol <-	selFDataCol[selFDataCol %in% names(fData(sqaPeptide$eset))] 
					
	cv <- sqaPeptide$cv
	names(cv) <- paste("cv",names(cv),sep="_")
	ratio <- sqaPeptide$ratio
	names(ratio) <- paste("log2ratio",names(ratio),sep="_")
	pValue <- sqaPeptide$pValue
	names(pValue) <- paste("pValue",names(pValue),sep="_")
	qValue <- sqaPeptide$qValue
	names(qValue) <- paste("qValue",names(qValue),sep="_")
	
	out <- cbind(
			fData(sqaPeptide$eset)[,selFDataCol]
			, exprs(sqaPeptide$eset)
			, cv
			, ratio	
			, pValue
			, qValue )[!fData(sqaPeptide$eset)$isFiltered,]
	
	write.table(out
			, file=userOptions$peptideTsvFilePath
			, sep="\t"
			, row.names=F
		)

	cat("INFO: CREATED FILE ", userOptions$peptideTsvFilePath,"\n")	
}	
	
if(exists("sqaProtein")){
	
	selFDataCol <- c("proteinName","proteinDescription","idScore","idQValue","nbPeptides","isNormAnchor")
	selFDataCol <- selFDataCol[selFDataCol %in%  names(fData(sqaProtein$eset))] 
	
	cv <- sqaProtein$cv
	names(cv) <- paste("cv",names(cv),sep="_")
	ratio <- sqaProtein$ratio
	names(ratio) <- paste("log2ratio",names(ratio),sep="_")
	pValue <- sqaProtein$pValue
	names(pValue) <- paste("pValue",names(pValue),sep="_")
	qValue <- sqaProtein$qValue
	names(qValue) <- paste("qValue",names(qValue),sep="_")
	
	out <- cbind(
			fData(sqaProtein$eset)[,selFDataCol]
			, exprs(sqaProtein$eset)
			, cv
			, ratio	
			, pValue
			, qValue )[!fData(sqaProtein$eset)$isFiltered,]
	
	# add median top3
	if(exists("esetTop3")){
		tmpOut <- getSignalPerCondition(esetTop3)
		names(tmpOut) <- paste("top3",names(tmpOut),sep="_")
		out <- cbind(out,tmpOut)
	}
	
	# add iBAQ	
	if(exists("esetIBAQ")){
		tmpOut <- exprs(esetIBAQ)
		names(tmpOut) <- paste("iBAQ",names(tmpOut),sep="_")
		out <- cbind(out,tmpOut)
	}
	
	write.table(out
			, file=userOptions$proteinTsvFilePath
			, sep="\t"
			, row.names=F
	)
	
	cat("INFO: CREATED FILE ", userOptions$proteinTsvFilePath,"\n")	
} 

### TSV EXPORT END

### EXPORT PARAMS
write.table(data.frame(
					param=row.names(data.frame(unlist(userOptions[names(userOptions)])))
					,value=as.vector((unlist(userOptions[names(userOptions)])))
			)
	,file=userOptions$paramsFilePath
	,sep="\t"
	,row.names=F
	,quote=F
)

cat("INFO: CREATED FILE ", userOptions$paramsFilePath,"\n")	

### EXPORT PARAMS

### EXPORT RDATA

if(userOptions$isSaveRObject){
	save.image(file=userOptions$rDataFilePath)
	cat("INFO: CREATED FILE ", userOptions$rDataFilePath,"\n")	
}
### EXPORT RDATA END


