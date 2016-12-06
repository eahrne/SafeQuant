#!/usr/local/bin/Rscript

###
# 0) Update SHEBANG ('#!/usr/bin/Rscript') to match location of your R installation
# 1) INSTALL SafeQuant
#	- or SET sqDirPath 
# 2) INSTALL PACKAGES
#	affy
#	limma
#	gplots
#	seqinr
#	corrplot
#	optparse
#	data.table

# Author: ahrnee-adm
###############################################################################

# TEST FILE
# /Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/new/peptides1_FILTERED.csv /Volumes/pcf01\$/Schmidt_Group/Databases/SwissProt_Databases/s_human_d_201405.fasta
# /Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/new/proteins1.csv
# /Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/new/TMT_6-Plex_Scaffold_Raw_Export_Example.xls

############################################################### INIT ############################################################### 
#### DEPENDANCIES
suppressWarnings(suppressPackageStartupMessages(library("affy", quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library("limma", quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(gplots, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(seqinr, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(corrplot, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(optparse, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(data.table, quiet=T)))

sourceDirOSX <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/R/"
sourceDirTPP <-  "/import/bc2/home/pcf/ahrnee/R/SafeQuant/R/"

# first check if dev or tpp mode
if(file.exists(sourceDirOSX) | file.exists(sourceDirTPP)){
	
	sourceDir <- ifelse(file.exists(sourceDirOSX),sourceDirOSX,sourceDirTPP)
	
	source(paste(sourceDir,"ExpressionAnalysis.R",sep=""))
	source(paste(sourceDir,"SafeQuantAnalysis.R",sep=""))
	source(paste(sourceDir,"Graphics.R",sep=""))
	source(paste(sourceDir,"IdentificationAnalysis.R",sep=""))
	source(paste(sourceDir,"Parser.R",sep=""))
	source(paste(sourceDir,"TMT.R",sep=""))
	source(paste(sourceDir,"UserOptions.R",sep=""))
	
}else if("SafeQuant" %in%  installed.packages()[,1]){ # used installed SafeQuant
	
	cat("Loading SafeQuant Library \n")
	library("SafeQuant")
	
}else{
	stop("SafeQuant Package not installed\n")
}

VERSION <- "2.3.1"

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

############################################################### PARSING ############################################################### 

if(userOptions$verbose) cat("PARSING INPUT FILE \n")	

# get file type
fileType <- .getFileType(userOptions$inputFile)

### Progenesis Export
if(fileType %in% c("ProgenesisProtein","ProgenesisFeature","ProgenesisPeptide")){
	
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
		
	}else if(fileType == "ProgenesisPeptide"){
		
		#"ProgenesisFeature"
		cat("INFO: PARSING PROGENESIS PEPTIDE EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisPeptideMeasurementCsv(file=userOptions$inputFile,expDesign=expDesign, uniqueProteins=userOptions$FUniquePeptides)
		
	}else{ 	#"ProgenesisFeature"
		
		cat("INFO: PARSING PROGENESIS FEATURE EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisFeatureCsv(file=userOptions$inputFile,expDesign=expDesign)
	}
	
# Scaffold Export (TMT data)
}else if(fileType == "ScaffoldTMT"){
	cat("INFO: PARSING SCAFFOLD RAW EXPORT FILE ",userOptions$inputFile, "\n" )
	
	# get default experimental design
	# six plex or ten plex ?
	# use default experimental design unless specified by the user
	nbPlex <- .getNbPlex(userOptions$inputFile) 
	if(nbPlex == 6){
		# 6-plex default: 1,2,3:4,5,6 
		expDesign <- data.frame(condition=paste("Condition",sort(rep(c(1,2),3)),sep=""),isControl=sort(rep(c(T,F),3),decreasing=T) )
	}else{
		# 10-plex default is "1,4,7,10:2,5,8:3,6,9"
		expDesign <- data.frame(condition=paste("Condition",c(1,2,3,1,2,3,1,2,3,1),sep=""),isControl=c(T,F,F,T,F,F,T,F,F,T) )
	}
	
	eset <- parseScaffoldRawFile(file=userOptions$inputFile,expDesign=expDesign)
	
	if(userOptions$TAdjustRatios){
		#if((nbPlex == 10)){ # only possible for tmt-10 plex
			nbCalMixSpectra <- sum( (fData(eset)$proteinName %in% names(CALIBMIXRATIOS)))
			if(nbCalMixSpectra < 100) stop("Not enough Calibration Mix spectra were found: ",nbCalMixSpectra, "\n ")
			cat("INFO: FOUND  ", nbCalMixSpectra ," Calibration Mix spectra\n")
			esetCalibMix <- .getCalibMixEset(eset)
			
			# discard calibration mix proteins
			eset <- eset[!(fData(eset)$proteinName %in% names(CALIBMIXRATIOS)),]
			
			intAdjObj <- .intensityAdjustment(eset, esetCalibMix)
			
		#}else{
		#	stop("Ratio Correction Not implemented for TMT 6-plex")
		#}
	}
	
	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)
	}
	
#	eset <- parseScaffoldRawFile(file=userOptions$inputFile,expDesign=expDesign)
	# apply specified experimental design
	eset <- createExpressionDataset(expressionMatrix=exprs(eset)[,rownames(expDesign)],expDesign=expDesign,featureAnnotations=fData(eset))
	
	if(!is.na(userOptions$scaffoldPTMSpectrumReportFile)){
		
		cat("INFO: ADDING SCAFFOLD PTM ANNOTATIONS \n")
		eset <- addScaffoldPTMFAnnotations(eset,userOptions$scaffoldPTMSpectrumReportFile)
		
	}
	
}else if(fileType == "MaxQuantProteinGroup"){
	
	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,data.frame(condition=paste("Condition",1:1000),isControl=c(F,1000)) )
	}else{
		stop("Please Specify Experimental Design")
	}
	eset <- parseMaxQuantProteinGroupTxt(userOptions$inputFile,expDesign=expDesign, method="auc")
	
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
	
	# dirty fix check if ACs in Progenesis file are stripped
	if(isStrippedACs(sample(fData(eset)$proteinName,100))){
		cat("INFO: RE-FORMATTING ACCESSION NUMBERS\n")
		names(proteinDB) <- stripACs(names(proteinDB))
	} 
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

# do not filter TMT data
if("pMassError" %in% names(fData(eset))  &&  (fileType != "ScaffoldTMT") ){
	### applicable to Progenesis feature Exports	
	
	if(is.na(userOptions$precursorMassFilter)){ # if not user specified  
		# automatically get precursor limits, X * sd of 50% top scoring
		userOptions$precursorMassFilter <- getMeanCenteredRange(fData(eset)$pMassError[fData(eset)$idScore > quantile(fData(eset)$idScore,na.rm=T)[3]],nbSd = 3)
		filter <- cbind(filter, pMassError=
						(fData(eset)$pMassError < userOptions$precursorMassFilter[1])
						| (fData(eset)$pMassError > userOptions$precursorMassFilter[2]) # precursor mass tolerance
		)
	}
}

if("ptm" %in% names(fData(eset))){
	
	# add motif-X and ptm coordinates
	if(exists("proteinDB")){
		
		cat("INFO: EXTRACTING PTM COORDINATES AND MOTIFS\n")
		#format 1) progensis  2) scaffold
		eset <- .addPTMCoord(eset,proteinDB,motifLength=6, isProgressBar=T,format= (fileType == "ScaffoldTMT") +1)
		
	}
	filter <- cbind(filter
			, ptm = !(grepl(userOptions$selectedModifName,as.character(fData(eset)$ptm),ignore.case=T))
			, nbPtmsPerPeptide = (fData(eset)$nbPtmsPerPeptide > userOptions$maxNbPtmsPerPeptide) )
	
}

if("peptide" %in% names(fData(eset))){
	filter <- cbind(filter
			, peptideLength =nchar(as.character(fData(eset)$peptide)) < userOptions$minPeptideLength 
			, charge =  fData(eset)$charge == 1 # discard singly charged
	)
}

if(!("nbPeptides" %in% names(fData(eset)))){
	### set nb peptides per protein
	eset <- setNbPeptidesPerProtein(eset)
}

filter <- cbind(filter,nbPeptides=(fData(eset)$nbPeptides < userOptions$minNbPeptidesPerProt))

# do not filter TMT data
#if(("idScore" %in% names(fData(eset))) && (fileType != "ScaffoldTMT")){
if(("idScore" %in% names(fData(eset)))){
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

# create paired experiemntal design()
if(userOptions$ECorrelatedSamples){
	eset <- createPairedExpDesign(eset)
}

# add filters etc to adjusted expressionSet
# update expDesign of intAdjObj$esetAd
if(exists("intAdjObj")){
	fData(intAdjObj$esetAdj) <- fData(eset)
	pData(intAdjObj$esetAdj) <- pData(eset)
	exprs(intAdjObj$esetAdj) <- exprs(intAdjObj$esetAdj)[,colnames(exprs(eset))]
} 


### non-pairwise stat test
statMethod <- c("")
if(userOptions$SNonPairWiseStatTest) statMethod <- c("all")
if(userOptions$SRawDataAnalysis){ # No Normalization
	esetNorm <- eset
	if(exists("intAdjObj")) intAdjObj$esetAdjNorm <- intAdjObj$esetAdj
}else{
	method <- c("global","median")
	# norm based on sum if norm anchor is specified 
	if(sum(fData(eset)$isNormAnchor) < nrow(eset)) method <- c("global","sum")
	esetNorm <- sqNormalize(eset, method=method)
	if(exists("intAdjObj")) intAdjObj$esetAdjNorm <- sqNormalize(intAdjObj$esetAdj, method=method )
}

### add pseudo (baseline) intensity
baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(esetNorm)[,1])),promille=5)
exprs(esetNorm)[  is.na(exprs(esetNorm)) | (exprs(esetNorm) <= 0)  ] <- 0
exprs(esetNorm) <- exprs(esetNorm) + baselineIntensity
if(exists("intAdjObj")){
	baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(intAdjObj$esetAdjNorm)[,1])),promille=5)
	exprs(intAdjObj$esetAdjNorm)[  is.na(exprs(intAdjObj$esetAdjNorm)) | (exprs(intAdjObj$esetAdjNorm) <= 0)  ] <- 0
	exprs(intAdjObj$esetAdjNorm) <- exprs(intAdjObj$esetAdjNorm) + baselineIntensity
}

if((fileType == "ProgenesisProtein") |  (fileType == "MaxQuantProteinGroup")){
	
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered  | isDecoy(fData(esetNorm)$proteinName)
	sqaProtein <- safeQuantAnalysis(esetNorm, method=statMethod, fcThrs=userOptions$ratioCutOff)
}else if((fileType == "ScaffoldTMT") && is.na(userOptions$scaffoldPTMSpectrumReportFile)){
	
	# roll-up protein level
	cat("INFO: ROLL-UP PROTEIN LEVEL\n")
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered |  isDecoy(fData(esetNorm)$proteinName)
	
	# correct TMT ratios
	if(userOptions$TAdjustRatios){
		fData(intAdjObj$esetAdjNorm)$isFiltered <- fData(esetNorm)$isFiltered
		intAdjObjProt <- intAdjObj
		intAdjObjProt$esetAdjNorm <- rollUp(intAdjObj$esetAdjNorm,featureDataColumnName= c("proteinName"))
		sqaProtein <- safeQuantAnalysis(rollUp(esetNorm,featureDataColumnName= c("proteinName")), method=statMethod,intensityAdjustmentObj=intAdjObjProt, fcThrs=userOptions$ratioCutOff )
	}else{
		sqaProtein <- safeQuantAnalysis(rollUp(esetNorm,featureDataColumnName= c("proteinName")), method=statMethod , fcThrs=userOptions$ratioCutOff)
	}
	
	fData(sqaProtein$eset)$isFiltered <- fData(sqaProtein$eset)$isFiltered | isDecoy(fData(sqaProtein$eset)$proteinName) | (fData(sqaProtein$eset)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	
}else{
	
	# roll-up peptide level
	cat("INFO: ROLL-UP PEPTIDE LEVEL\n")
	
	# correct TMT ratios
	if(userOptions$TAdjustRatios){
		cat("WARN: Ratio correction not yet implemented in this anlysis mode \n")			
	}
	
	esetPeptide <- rollUp(esetNorm,featureDataColumnName= c("peptide","ptm"))
	
	# fdr filter
	# replace qValues by rollUp level qValues ()
	esetPeptide <- addIdQvalues(esetPeptide)
	
	if(fileType == "ScaffoldTMT"){
		fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | (fData(esetPeptide)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	}else{
		# update filter to exclude peptide level hight qValues
		fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | (fData(esetPeptide)$idQValue > userOptions$fdrCutoff) | (fData(esetPeptide)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	}
	
	if(userOptions$proteinQuant){
		cat("INFO: ROLL-UP PROTEIN LEVEL\n")
		esetProtein <- rollUp(esetPeptide,featureDataColumnName= c("proteinName"))
		esetProtein <- addIdQvalues(esetProtein)
		
		if(fileType == "ScaffoldTMT"){
			fData(esetProtein)$isFiltered <- fData(esetProtein)$isFiltered | isDecoy(fData(esetProtein)$proteinName) | (fData(esetProtein)$nbPeptides <  userOptions$minNbPeptidesPerProt)
		}else{
			fData(esetProtein)$isFiltered <- fData(esetProtein)$isFiltered | (fData(esetProtein)$idQValue > userOptions$fdrCutoff) | isDecoy(fData(esetProtein)$proteinName) | (fData(esetProtein)$nbPeptides <  userOptions$minNbPeptidesPerProt)
		}
		sqaProtein <- safeQuantAnalysis(esetProtein, method=statMethod, fcThrs=userOptions$ratioCutOff)
		
	}
	
	fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | isDecoy(fData(esetPeptide)$proteinName)
	sqaPeptide <- safeQuantAnalysis(esetPeptide, method=statMethod, fcThrs=userOptions$ratioCutOff)
	
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered | isDecoy(fData(esetNorm)$proteinName) | (fData(esetNorm)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	
	if(userOptions$top3 & userOptions$proteinQuant){
		cat("INFO: ROLL-UP TOP3\n")
		esetTop3 <-  rollUp(esetPeptide,featureDataColumnName= c("proteinName"), method="top3")
	}
}

### IBAQ
if(userOptions$iBAQ & userOptions$proteinQuant){
	cat("INFO: CALCULATING IBAQ VALUES\n")
	if(exists("proteinDB")){
		esetIBAQ <-  getIBAQEset(sqaProtein$eset, proteinDB=proteinDB)
	}else{
		cat("ERROR: proteinDB NOT FOUND NO iBAQ VALUES CALCULATED\n")
	}
}

### EXPRESSION ANALYSIS END

############################################################### EXPORTS ############################################################### 
cat("INFO: PREPARING EXPORTS","\n")

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

############################### GRAPHICS

# plot protein or peptide level results
if(exists("sqaProtein")){
	sqaDisp <- sqaProtein
	lab <- "Protein"
}else{
	sqaDisp <- sqaPeptide
	lab <- "Peptide"
}
### only disp. a subset for some plots
rowSelEset <- 1:nrow(eset) %in% sample(nrow(eset),min(c(2000,nrow(eset))) ,replace=F)
rowSelSqaDisp <- 1:nrow(sqaDisp$eset) %in% sample(nrow(sqaDisp$eset),min(c(2000,nrow(sqaDisp$eset))) ,replace=F)

pdf(userOptions$pdfFile)
parDefault <- par()
CONDITIONCOLORS <- .getConditionColors(esetNorm)

### EXPDESIGN PLOT
plotExpDesign(esetNorm, version=VERSION)
### EXPDESIGN PLOT END

### IDENTIFICATION PLOTS
if(userOptions$verbose) cat("INFO: IDENTIFICATION PLOTS \n")
par(mfrow=c(2,2))

#.idOverviewPlots()
#@ NOT CRAN COMPATIBLE	
.idOverviewPlots(userOptions=userOptions
		,esetNorm=esetNorm
		,fileType=fileType
		,sqaPeptide= ifelse(exists("sqaPeptide"),list(sqaPeptide),list(NA))[[1]]# HACK to pass check
		,sqaProtein= ifelse(exists("sqaProtein"),list(sqaProtein),list(NA))[[1]] # HACK to pass check
)

if(fileType %in% c("ProgenesisFeature","ProgenesisPeptide")){
	par(mfrow=c(3,2))
	.idPlots(eset, selection=c(1,3), main="Feature Level", qvalueThrs=userOptions$fdrCutoff, userOptions=userOptions)
	if(exists("sqaPeptide")) .idPlots(sqaPeptide$eset, selection=c(1,3), main="Peptide Level", qvalueThrs=userOptions$fdrCutoff)
	if(exists("sqaProtein")) .idPlots(sqaProtein$eset, selection=c(1,3), main="Protein Level", qvalueThrs=userOptions$fdrCutoff)
}	
par(parDefault)
### IDENTIFICATIONS PLOTS END
### QUANT. QC PLOTS 
if(userOptions$verbose) cat("INFO: QUANT QC. PLOTS \n")

### MASS ERROR
par(parDefault)
#if("pMassError" %in% names(fData(eset))){
if(fileType %in% c("ProgenesisFeature","ProgenesisPeptide")){	
	par(mfrow=c(2,1), mar=c(4.5,6.1,4.1,6.1))
	plotPrecMassErrorDistrib(eset, pMassTolWindow=userOptions$precursorMassFilter)
	
	plotPrecMassErrorVsScore(eset[rowSelEset,], pMassTolWindow=userOptions$precursorMassFilter)
	par(parDefault)
}


layout(rbind(c(1,2), c(3,3)))
### missing values
missinValueBarplot(eset)

### total intensity sum
barplotMSSignal(eset)
par( mar=c(6.5,5.1,2.5,3.1))
cvBoxplot(sqaDisp$eset)

par(parDefault)

### CORRELATION PLOTS
### COR OR PAIRS PLOT. IF FEWER THAN X SAMPLES 

if(ncol(sqaDisp$eset) < 8){
	pairsAnnot(log10(exprs(sqaDisp$eset))[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered ,],textCol=as.character(CONDITIONCOLORS[pData(sqaDisp$eset)$condition,]))
}else{
	.correlationPlot(log10(exprs(sqaDisp$eset))[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,], labels=as.character(unique(pData(sqaDisp$eset)$condition)), textCol=as.character(CONDITIONCOLORS[pData(sqaDisp$eset)$condition,]))
}
### COR OR PAIRS PLOT. IF FEWER THAN X CONDITIONS 
if(length(unique(pData(sqaDisp$eset)$condition)) < 8){
	pairsAnnot(log10(getSignalPerCondition(sqaDisp$eset[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,]))[,as.character(unique(pData(sqaDisp$eset)$condition)) ],textCol=as.character(CONDITIONCOLORS[as.character(unique(pData(sqaDisp$eset)$condition)),]))
}else{
	.correlationPlot(log10(getSignalPerCondition(sqaDisp$eset[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,]))[,as.character(unique(pData(sqaDisp$eset)$condition)) ],textCol=as.character(CONDITIONCOLORS[as.character(unique(pData(sqaDisp$eset)$condition)),]))
}

par(parDefault)


### TMT calibration mix
if(exists("intAdjObj")){ 
	
	# plot adjusted ratios vs org ratio
	# boxplot noise fraction
	#if(ncol(sqaDisp$ratio) > 1) par(mfrow=c(2,2))	
	boxplot(intAdjObj$noiseFraction*100, border=ifelse(intAdjObj$selectedPairs,"blue","black")
			, ylab="Noise Fraction (%)",xlab="Calibration Mix Pair", cex.axis=1.5,cex.lab=1.5)
	
	if(ncol(sqaDisp$ratio) > 1) par(mfrow=c(2,2))	
	plotAdjustedVsNonAdjustedRatio(sqaDisp$ratio,sqaDisp$unAdjustedRatio)
	par(parDefault)

}



### QUANT. QC PLOTS END

par(parDefault)
if(userOptions$verbose) cat("INFO: HEAT MAP \n")
hClustHeatMap(sqaDisp$eset[(1:nrow(sqaDisp$eset) %in% sample(nrow(sqaDisp$eset),min(c(nrow(sqaDisp$eset),10000)))) &  !fData(sqaDisp$eset)$isFiltered,],main= paste(lab,"Level"))

### QUANT. STAT. PLOTS 

### VAILD FEATURES VS. pValue/qValue
if(userOptions$verbose) cat("INFO: QUANT RES. PLOTS \n")

par(mfrow=c(1,2))
plotNbValidDeFeaturesPerFDR(sqaDisp,
		upRegulated=F
		,log2RatioCufOff=log2(userOptions$ratioCutOff)
		,pvalRange=c(0,0.15)
		,pvalCutOff=userOptions$deFdrCutoff
		,isLegend=T
		,isAdjusted=T
		,ylab=paste(lab, "Counts")
		,main="DOWN REGULATION"
)

plotNbValidDeFeaturesPerFDR(sqaDisp,
		upRegulated=T
		,log2RatioCufOff=log2(userOptions$ratioCutOff)
		,pvalRange=c(0,0.15)
		,pvalCutOff=userOptions$deFdrCutoff
		,isLegend=F
		,isAdjusted=T
		,ylab=paste(lab, "Counts")
		,main="UP REGULATION"
)

par(parDefault)

plotVolcano(sqaDisp
		, main=paste(lab,"Level")
		, ratioThrs= userOptions$ratioCutOff
		, pValueThreshold= userOptions$deFdrCutoff
		, adjusted = T)


if(userOptions$eBayes){
	
	par(mfrow=c(1,2))
	plotNbValidDeFeaturesPerFDR(sqaDisp,
			upRegulated=F
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=T
			,isAdjusted=F
			,ylab=paste(lab, "Counts")
			,main="DOWN REGULATION"
	)
	
	plotNbValidDeFeaturesPerFDR(sqaDisp,
			upRegulated=T
			,log2RatioCufOff=log2(userOptions$ratioCutOff)
			,pvalRange=c(0,0.15)
			,pvalCutOff=userOptions$deFdrCutoff
			,isLegend=F
			,isAdjusted=F
			,ylab=paste(lab, "Counts")
			,main="UP REGULATION"
	)
	par(parDefault)
	
	plotVolcano(sqaDisp
			, main=paste(lab,"Level")
			, ratioThrs= userOptions$ratioCutOff
			, pValueThreshold= userOptions$deFdrCutoff
			, adjusted = F)
	
	par(mfrow=c(2,2))
	if(nrow(CONDITIONCOLORS) > 4) par(mfrow=c(3,3))
	.allpValueHist(sqaDisp)
	plotQValueVsPValue(sqaDisp, lim=c(0,1))
	par(parDefault)
}		

### QUANT. STAT. PLOTS END

par(parDefault)


### SOME ADDITIONAL QC PLOTS

if(userOptions$addQC){
	
#	if(exists("sqaPeptide")){
#		plotXYDensity(fData(sqaPeptide$eset)$retentionTime,fData(sqaPeptide$eset)$pMassError, disp=c("")
#			, xlab="Retention time (min)"
#			, ylab="Precursor Mass Error (ppm)"
#			, cex.axis=1.5
#			, cex.lab=1.5)
	
	if( all(c("retentionTime","pMassError")  %in% names(fData(eset)) )){
		plotXYDensity(fData(eset)$retentionTime,fData(eset)$pMassError, disp=c("")
			, xlab="Retention time (min)"
			, ylab="Precursor Mass Error (ppm)"
			, cex.axis=1.5
			, cex.lab=1.5)
	
		abline(h=c(userOptions$precursorMassFilter[1],0,userOptions$precursorMassFilter[2]),lty=2, lwd=2)
		
		# rt vs signal
		sel <- 1:nrow(esetNorm) %in% sample(nrow(esetNorm),min(c(4000,nrow(esetNorm))) ,replace=F) & (!(fData(esetNorm)$proteinName %in% names(CALIBMIXRATIOS)))
		plotRTNormSummary(esetNorm[sel,])
		
		par(mfrow=c(2,2))
		plotRTNorm(getRTNormFactors(esetNorm[sel,], minFeaturesPerBin=100),esetNorm[sel,])
	
		par(parDefault)
		
	}
	
	par(mfrow=c(2,2))
	#all ma plots
	for(s in colnames(exprs(esetNorm))){
		sel <- 1:nrow(esetNorm) %in% sample(nrow(esetNorm),min(c(4000,nrow(esetNorm))) ,replace=F) & (!(fData(esetNorm)$proteinName %in% names(CALIBMIXRATIOS)))
		maPlot(esetNorm[sel,],sample=s)
	}
	
	
	
	par(parDefault)
	
}



cat("INFO: CREATED FILE ", userOptions$pdfFile,"\n")

graphics.off()

############################### GRAPHICS END

### TSV EXPORT

if(exists("sqaPeptide")){ 
	
	selFDataCol <- c("peptide","proteinName","proteinDescription", "idScore","idQValue"
			,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures" ) 
	selFDataCol <-	selFDataCol[selFDataCol %in% names(fData(sqaPeptide$eset))] 
	
	### add modif coord
	if("motifX" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"motifX","modifCoord")
	}
	
	### add allAccessions
	if("allAccessions" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"allAccessions")
	}
	
	### add ptmPeptide
	if("ptmPeptide" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmPeptide")
	}
	
	### add ptmLocProb
	if("ptmLocProb" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmLocProb")
	}
	
	### add ptmLocMascotConfidence
	if("ptmLocMascotConfidence" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmLocMascotConfidence")
	}
	
	cv <- sqaPeptide$cv
	names(cv) <- paste("cv",names(cv),sep="_")
	ratio <- sqaPeptide$ratio
	
	if(length(names(ratio)) > 0 ) names(ratio) <- paste("log2ratio",names(ratio),sep="_")
	
	pValue <- sqaPeptide$pValue
	if(length(names(pValue)) > 0 ) names(pValue) <- paste("pValue",names(pValue),sep="_")
	qValue <- sqaPeptide$qValue
	if(length(names(qValue)) > 0 )  names(qValue) <- paste("qValue",names(qValue),sep="_")
	
	medianSignalDf <- getSignalPerCondition(sqaPeptide$eset)
	names(medianSignalDf) <- paste("medianInt",names(medianSignalDf),sep="_")
	
	out <- cbind(
			fData(sqaPeptide$eset)[,selFDataCol]
			, exprs(sqaPeptide$eset)
			, medianSignalDf
			, cv
			, ratio	
			, pValue
			, qValue )[!fData(sqaPeptide$eset)$isFiltered,]
	
	### add unadjusted ratios if TMT ratio correction
	if(userOptions$TAdjustRatios){
		unadjPeptideRatios <- sqaPeptide$unAdjustedRatio[!fData(sqaPeptide$eset)$isFiltered,]
		names(unadjPeptideRatios) <- paste("log2_unadjRatio",names(sqaPeptide$ratio),sep="_")
		out <- cbind(out,unadjPeptideRatios)
	}
	
	### paired expDesign ratio export
	if("subject" %in% names(pData(sqaPeptide$eset))){
		allRatios <- getRatios(sqaPeptide$eset,method="paired")[!fData(sqaPeptide$eset)$isFiltered,]
		names(allRatios) <- paste("log2_pairedRatio",names(allRatios),sep="_")
		out <- cbind(out,allRatios)
	}
	
	write.table(out
			, file=userOptions$peptideTsvFilePath
			, sep="\t"
			, row.names=F
	)
	
	cat("INFO: CREATED FILE ", userOptions$peptideTsvFilePath,"\n")	
}	

if(exists("sqaProtein")){
	
	selFDataCol <- c("proteinName","proteinDescription","idScore","idQValue","nbPeptides")
	selFDataCol <- selFDataCol[selFDataCol %in%  names(fData(sqaProtein$eset))] 
	
	### add allAccessions
	if("allAccessions" %in% names(fData(sqaProtein$eset))){
		selFDataCol <- c(selFDataCol,"allAccessions")
	}
	
	cv <- sqaProtein$cv
	names(cv) <- paste("cv",names(cv),sep="_")
	ratio <- sqaProtein$ratio
	if(length(names(ratio)) > 0 ) names(ratio) <- paste("log2ratio",names(ratio),sep="_")
	pValue <- sqaProtein$pValue
	if(length(names(pValue)) > 0 ) names(pValue) <- paste("pValue",names(pValue),sep="_")
	qValue <- sqaProtein$qValue
	if(length(names(qValue)) > 0 ) names(qValue) <- paste("qValue",names(qValue),sep="_")
	
	medianSignalDf <- getSignalPerCondition(sqaProtein$eset)
	names(medianSignalDf) <- paste("medianInt",names(medianSignalDf),sep="_")
	
	out <- cbind(
			fData(sqaProtein$eset)[,selFDataCol]
			, exprs(sqaProtein$eset)
			, medianSignalDf
			, cv
			, ratio	
			, pValue
			, qValue )[!fData(sqaProtein$eset)$isFiltered,]
	
	# add median top3
	if(exists("esetTop3")){
		
		# medians
#		tmpOut <- getSignalPerCondition(esetTop3)
#		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
#		names(tmpOut) <- paste("medianInt_top3",names(tmpOut),sep="_")
		
		tmpOut <- exprs(esetTop3)
		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
		colnames(tmpOut) <- paste("top3",colnames(tmpOut),sep="_")
				
		out <- cbind(out,tmpOut)
	}
	
	# add iBAQ	
	if(exists("esetIBAQ")){
		
		# medians 
#		tmpOut <- getSignalPerCondition(esetIBAQ)
#		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
#		names(tmpOut) <- paste("medianInt_top3",names(tmpOut),sep="_")
		
		tmpOut <- exprs(esetIBAQ)
		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
		colnames(tmpOut) <- paste("iBAQ",colnames(tmpOut),sep="_")
		
		out <- cbind(out,tmpOut)
		
	}
	
	### add unadjusted ratios if TMT ratio correction
	if(userOptions$TAdjustRatios){
		unadjProteinRatios <- sqaProtein$unAdjustedRatio[!fData(sqaProtein$eset)$isFiltered,]
		names(unadjProteinRatios) <- paste("log2_unadjRatio",names(sqaProtein$ratio),sep="_")
		out <- cbind(out,unadjProteinRatios)
	}
	
	
	### paired expDesign ratio export
	if("subject" %in% names(pData(sqaProtein$eset))){
		allRatios <- getRatios(sqaProtein$eset,method="paired")[!fData(sqaProtein$eset)$isFiltered,]
		names(allRatios) <- paste("log2_pairedRatio",names(allRatios),sep="_")
		out <- cbind(out,allRatios)
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
