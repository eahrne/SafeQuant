#!/usr/bin/Rscript

# TODO: Add comment
# 
# Author: erikahrne
###############################################################################
#### LOAD FUNCTIONS, VERSION 


### get directory of script
initial.options <- commandArgs(trailingOnly = FALSE)
script.basename  <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
source(paste(script.basename,"/functions/ProgenesisPPFunctions.R",sep=""))
print(paste(script.basename,"/config/VERSION.R",sep=""))
source(paste(script.basename,"/config/VERSION.R",sep=""))

###

### CMD OPTIONS
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
				help="Print extra output [default %default]"),
		make_option(c("-i", "--inputFile"), type="character", default="",
				help="Progenesis csv input file (REQUIRED)",
		),	
		make_option(c("-o", "--outputDir"), type="character", default="",
				help="Results Output Directory [default ./]",
		),
		make_option(c("-l", "--resultsFileLabel"), type="character", default="analysisResults",
				help=" .pdf & .csv file labels (prefix) [default %default]", 
		),
		make_option(c("-a", "--proteinAccesionFilter"), type="character", default=".",
				help="Filter features by Accession [default %default] (all features kept)",
				metavar="protein accession regexpr"),
		make_option(c("-t", "--modificationTypeFilter"), type="character", default=".",
				help="Filter features by Modification Type [default %default] (all features kept)",
				metavar="modification name regexpr"),
#		make_option(c("-s", "--scoreMassDiffFilter"), action="store_true", default=FALSE,
#				help="Filter out feature with (both) score < YY & massDiff > XX  [default %default]"),
		
		make_option(c("-m", "--precursorMassFilter"), type="double", default=10,
				help="Precursor mass filter 
						min 1 ppm [default %default ppm]"),
		
		make_option(c("-c", "--controlCondition"), type="integer", default=1,
				help="Specifiy condition number to be used as control [default %default]",
				metavar="selected control condition"),
						
		make_option(c("-f", "--fdrCutoff"), type="double", default=0.01,
				help="Calculating intetnsity ratios for peptides with 
	a score corresponding to a fdr cut-off < specified value. [0-1] [default %default]",
				metavar="Peptide level FDR cutoff"),
		make_option(c("-r", "--ratioCutoff"), type="double", default=1.5,
				help="Intensity ratio (fold change) cut-off 
	used for graphics and results export. >1 [default %default]",
				metavar="Intensity ratio cutoff"),
		make_option(c("-b", "--logOddsCutoff"), type="double", default=2.2,
				help="log odds cut-off used for graphics.
	High-lighting features with a log odds ratio > specified value. [default %default]",
				metavar="Intensity ratio logOdds cutoff"),
		make_option(c("-q", "--deFdrCutoff"), type="double", default=0.01,
				help="fdr cut-off used for graphics.
	High lighting features with a qval < specified value. [0-1] [default %default]",
				metavar="Diff. expr. fdr cut-off"),
#		make_option(c("-m", "--meanFC"), action="store_true", default=FALSE,
#				help="Display ratio of mean intensities per condition instead of median based intensity ratio [default %default]"),
		make_option(c("-x", "--topX"), type="integer", default=0,
				help="creates csv output including the top X 
	most intense features per protein [default %default = off]",
				metavar="top X peptides per protein"),
		make_option(c("-n", "--naReplace"), type="double", default=100,
				help="Intensity replacing missing values. [default %default]",
		),	
		make_option(c("-d", "--proteaseTarget"), type="character", default="KR",
				help=" protease target aa's [default (Trypsin) %default]", 
		),
		make_option(c("-y", "--test"), action="store_true", default=FALSE,
				help="TEST (parse only first 1000 entries)  [default %default]"),
		
		make_option(c("-s", "--saveRObject"), action="store_true", default=FALSE,
				help="Save R objects in 'label'.RData file [default %default]"),
		
		make_option(c("-z", "--selectedGraphics"), type="character", default="",
				help="Excluded Graphics: give letter for each plot to exclude ex: -z iv 
		(creates all plots but intensity density plots & volcano plot)
			experimental design (include e)
			protein score distribution related plots (include f)
			intensity density plots (include i)
			
			volcano plot(s) (include v)
			sample and condition intensity correlation plots (include c)
			sample MVA plots (include m)
			hierarchical clustering plots (include h)	
			differential expression fdr plot (include d)	
			qval based volcano plot instead of log odds based (include l)	
			[default (all plots) %default]")

)

#spectral counts density plot (include s)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser( prog=paste("SafeQuant",VERSION), option_list=option_list))

#progenesisFile <- '/Users/erikahrne/Documents/Projects/progenesis_postProcess/proteins2.csv'
#progenesisFile <- '/Users/erikahrne/Documents/Projects/progenesis_postProcess/tmp/proteins.csv'
progenesisFile <- opt$inputFile
if( progenesisFile == "" | !file.exists(progenesisFile)){
	print(paste("ERROR. Please specify input file.",progenesisFile, "Not found!"))
	q(status=-1)
}

outputDir <- opt$outputDir
if(!file.exists(outputDir) & outputDir != "" ){
	print(paste("ERROR. No such directory",outputDir))
	q(status=-1)
}else if(substr(outputDir,nchar(outputDir),nchar(outputDir)) != "/"){
	if(opt$verbose){
		print(paste("added slash to outputdir")) 
	}
	if(!(outputDir== "")){
		outputDir <- paste(outputDir,"/", sep="")
	}
}

resultsFileLabel <- opt$resultsFileLabel

#pdfFile <- '/Users/erikahrne/dev/R/workspace/ProgenesisPP/out/bioC/bioctest.pdf'
#csvFile <- '/Users/erikahrne/dev/R/workspace/ProgenesisPP/out/bioC/bioctest.csv'
pdfFile <- paste(outputDir,resultsFileLabel,".pdf",sep="")
csvFile <- paste(outputDir,resultsFileLabel,".csv",sep="")

selectedProteinName <- opt$proteinAccesionFilter
selectedModifName <- opt$modificationTypeFilter

scoreMassDiffFilter <- opt$scoreMassDiffFilter

precursorMassFilter <- opt$precursorMassFilter
if(is.na(opt$precursorMassFilter) | precursorMassFilter <= 0 ){
	print(paste("ERROR. precursorMassFilter must be > 1 ppm. You specified",precursorMassFilter)) 
	q(status=-1)
}

selectedControlCond <- opt$controlCondition
if(is.na(opt$controlCondition) | selectedControlCond < 1 ){
	print(paste("ERROR. controlCondition must be >= 1. You specified ", selectedControlCond))
	q(status=-1)
}

fdrCutoff <- opt$fdrCutoff
if(is.na(opt$fdrCutoff) | fdrCutoff <= 0 | fdrCutoff > 1 ){
	print(paste("ERROR. fdrCutOff must be in the range [0-1]. You specified",fdrCutoff)) 
	q(status=-1)
}

ratioCutOff <- opt$ratioCutoff
if(is.na(opt$ratioCutoff) | ratioCutOff <= 1){
	print(paste("ERROR. ratioCutoff must be > 1. You specified",ratioCutOff)) 
	q(status=-1)
}

BCutoff <- opt$logOddsCutoff
if(is.na(opt$logOddsCutoff)){
	print(paste("ERROR. Invalid log Odds Cutoff. You specified",BCutoff)) 
	q(status=-1)
}

deFdrCutoff <- opt$deFdrCutoff
if(is.na(opt$deFdrCutoff) | deFdrCutoff <= 0 | deFdrCutoff > 1 ){
	print(paste("ERROR. deFdrCutoff must be in the range [0-1]. You specified",deFdrCutoff)) 
	q(status=-1)
}

#meanFC = opt$meanFC

topX <- opt$topX
if(is.na(opt$topX) | topX < 0){
	print(paste("ERROR. topX must be >= 0. You specified",topX)) 
	q(status=-1)
}

topXCsvFile <- paste(outputDir,resultsFileLabel,"_top",topX,".csv",sep="")

minIntensity <- opt$naReplace
if(is.na(opt$naReplace) | minIntensity <= 0){
	print(paste("ERROR. deFdrCutoff must be > 0. You specified",minIntensity)) 
	q(status=-1)
}

proteaseTarget <- opt$proteaseTarget

test <- opt$test

############ grapihcs 

dispExpDesign <- !regexpr("e",opt$selectedGraphics) > -1
fdrPlots <- !regexpr("f",opt$selectedGraphics) > -1
intensityDensityPlots <- !regexpr("i",opt$selectedGraphics) > -1
#spectralcountsDensityPlots <- !regexpr("s",opt$selectedGraphics) > -1
volcanoPlots <- !regexpr("v",opt$selectedGraphics) > -1
mvaPlots <- !regexpr("m",opt$selectedGraphics) > -1
correlationPlots <- !regexpr("c",opt$selectedGraphics) > -1
hClustPlot <- !regexpr("h",opt$selectedGraphics) > -1
deFDRPlot <- !regexpr("d",opt$selectedGraphics) > -1
qvalVolcano <- regexpr("l",opt$selectedGraphics) > -1
#lodsDistributionPlots <- regexpr("l",opt$selectedGraphics) > -1

############ grapihcs end 

### CMD OPTIONS END

### LOAD SCRIPTS

#suppressPackageStartupMessages( library(MASS))
suppressPackageStartupMessages(library(limma))
#suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(gplots))
#source("/Users/erikahrne/dev/R/workspace/ProgenesisPP/functions/ProgenesisPPFunctions.R")



### LOAD SCRIPTS END

### PARSE CSV FILE

if(opt$verbose){
	print(paste("PARSING CSV FILE",progenesisFile))
}

	### READ EXPERIMENTAL DESIGN

line2header <- read.csv(file=progenesisFile,head=TRUE,sep=",",check.names=FALSE, skip=1,nrows=1, strip.white=TRUE)
line2header <- gsub(" ","",names(line2header))
line2header.names.unique <- unique(line2header)[1:(length(unique(line2header))-1)]
condition.names <- line2header.names.unique[nchar(line2header.names.unique) > 0]

### replace condition names with 1,2 ...
i <- 1
for(cond in condition.names){
	
	line2header[line2header == cond] <- as.character(i)
	i <- i+1
}
orgCondNames <- condition.names
condition.names <- as.character(1:length(condition.names))

expDesign <- list()
i <- 1
for(i in 1:length(line2header)){
	lab <- line2header[i]
	
	if(lab %in% condition.names){
		j <- i+1
		while(nchar(line2header[j]) == 0 ){
			j <- j+1
		}
		replicates <- j-i
		if(!(lab %in% names(expDesign))){
			expDesign[[lab]] <- replicates
		}
	}
}

nbSamples <- 0
for(cond in condition.names){
	nbSamples <- nbSamples+ expDesign[[cond]]
}

enoughReplicates <- sum(unlist(expDesign)) >= 2*length(names(expDesign))

if(opt$verbose){
	print("Experimental design")
	print(line2header)
	print(expDesign)
}

if(opt$verbose){
	print("condition.names")
	print(condition.names)
}

nbConditions <- length(condition.names)

	### READ EXPERIMENTAL DESIGN END

 ### avoid parsing the whole file


#if(opt$verbose){
	print("Reading all peptide entries.. this can take a while")
#}

if(test){
	data <- read.csv(file=progenesisFile,head=TRUE,sep=",", skip=2, nrows=13000)
}else{
	data <- read.csv(file=progenesisFile,head=TRUE,sep=",", skip=2)
}

data <- data[!is.na(data$Score),]

col.names <- names(data)[regexpr("X",names(data)) == -1]
	
if(opt$verbose){
	print("col.names")
	print(col.names)
	
	print("nb peptides")
	print(length(rownames(data)))
	
}

rawIntColumnStart.tmp <- 13 +nbSamples
rawIntColumnEnd.tmp <- rawIntColumnStart.tmp+nbSamples-1
spectralCountColumnsStart.tmp <-  rawIntColumnEnd.tmp +(2*nbSamples)+ 9	
spectralCountColumnsEnd.tmp <- spectralCountColumnsStart.tmp+nbSamples-1

### PARSE CSV FILE END

### CREATE DATA FRAME OF UNIQUE PEPTIDES (max score & min mass Error (ppm) selected )

petideIds <-  paste(data$Sequence,data$Variable.modifications...position..description. ,sep="_:")
peptidesIds.unique <- unique(petideIds)

peptides <- c()
modifs <- c()
proteins <- c()
descriptions =  c()
minMassDiffs <- c()
maxScores <-c()
spc <- c()
rawIntensities <- c()
summedEntries <- c()

#if(opt$verbose){
	print("Grouping features by unique peptide sequence.. this can also take a while")
#}

pbSum <- txtProgressBar(min = 0, max = length(peptidesIds.unique), style = 3)

i <-1
for(pept in peptidesIds.unique){
	
#	if(progressBar){
		setTxtProgressBar(pbSum, i)
		i <-  i+1
#	}
		
	indices = which(pept == petideIds )
	index1 = indices[1]
	
	peptides = c(peptides, as.character(data$Sequence[index1]))
	modifs = c(modifs,as.character(data$Variable.modifications...position..description.[index1]))
	proteins = c(proteins,as.character(data$Protein[index1]))
	descriptions =  c(descriptions,as.character(data$Description[index1]))
	minMassDiffs = c(minMassDiffs, min(data$Mass.error..ppm.[indices]))
	maxScores = c(maxScores,max(data$Score[indices]))
	summedEntries = c(summedEntries, length(indices))
	
	###sum up intensities
	intensities.tmp = data[indices,rawIntColumnStart.tmp:rawIntColumnEnd.tmp]
	sumRawIntensities = list()
	sumRawIntensities[[pept]] = apply(intensities.tmp,2,FUN=sum)
	rawIntensities = c(rawIntensities,sumRawIntensities) 
	
	###sum up spectral counts
	spC.tmp = data[indices,spectralCountColumnsStart.tmp:spectralCountColumnsEnd.tmp]
	sumSpC = list() 
	sumSpC[[pept]] = apply(spC.tmp,2,FUN=sum)
	spc = c(spc,sumSpC)
	
}
close(pbSum)
cat("\n")

### COUNT MIS-CLEAVAGE SITES

proteaseRegexpr <- paste("[",proteaseTarget,"]",".",sep="")
missedCleavageSites <- unlist(lapply(peptides,FUN= function(f){return(sum(unlist(gregexpr(proteaseRegexpr,f))>0))}))
### COUNT MIS-CLEAVAGE SITES END

### @TODO allow user to specify control condition
### reorder rawData, expDesign, orgCondNames, spectralCounts,  placing control condition in first position

rawData.tmp <- t(data.frame(rawIntensities))
spectralCounts.tmp <- t(data.frame(spc))  

if(selectedControlCond != 1){
	
	if(selectedControlCond <= nbConditions){
		rD.tmp <-data.frame( row.names= rownames(rawData.tmp))
		rD.control.tmp <- data.frame( row.names= rownames(rawData.tmp))
		
		expDesign.tmp <-list()
		
		orgCondNames.tmp <- c()
		
		sC.tmp <- data.frame( row.names= rownames(rawData.tmp))
		sC.control.tmp <- data.frame( row.names= rownames(rawData.tmp))
		
		
		i <-1
		j <-2
		for(cond in names(expDesign)){
			
			condNbSamples <-  expDesign[[cond]]
			cond <- as.numeric(cond)
			
			if(cond == selectedControlCond){
				rD.control.tmp <- rawData.tmp[,i:(i+condNbSamples-1)]
				sC.control.tmp <- spectralCounts.tmp[,i:(i+condNbSamples-1)]
				expDesign.tmp[["1"]] <- condNbSamples
				
			}else{
				rD.tmp <-data.frame(rD.tmp, rawData.tmp[,i:(i+condNbSamples-1)] )
				sC.tmp <-data.frame(sC.tmp, spectralCounts.tmp[,i:(i+condNbSamples-1)] )
				expDesign.tmp[[as.character(j)]] <- condNbSamples
				orgCondNames.tmp <- c(orgCondNames.tmp, orgCondNames[cond])
				
				
				j <- j+1
			}  
			
			i <-i+condNbSamples
			
		}
		
		### put control cond first 
		rD.tmp <- data.frame(rD.control.tmp,rD.tmp)
		sC.tmp <- data.frame(sC.control.tmp,sC.tmp)
		orgCondNames.tmp <- c(orgCondNames[selectedControlCond],orgCondNames.tmp)
		
		#sC.tmp <-  data.frame(spectralCounts.tmp[,selectedControlCond], sC.tmp)
		#names(sC.tmp) <- orgCondNames.tmp
#	print("NEW")
#	print(names(rD.tmp)) 
#	print(names(expDesign))
#	print(names(expDesign.tmp))
#	print(expDesign.tmp)
#	print(orgCondNames.tmp)
#	print(colnames(spectralCounts.tmp))
#	print(colnames(sC.tmp))
		
		### apply changes
		expDesign <- expDesign.tmp
		rawData.tmp <-  rD.tmp
		spectralCounts.tmp <- sC.tmp 
		orgCondNames <- orgCondNames.tmp
		
		
	}else{
		
		print(paste("ERROR. Invalid controlCondition must be <=", nbConditions , ", You specified", selectedControlCond))
		q(status=-1)
		
	}
}


dataSummed <- data.frame(peptides
			,modifs
			,proteins
			,descriptions		
			,minMassDiffs
			,maxScores
			,missedCleavageSites		
			,summedEntries
			,rawData.tmp
			,spectralCounts.tmp
			)

	
rawIntColumnStart <- 9
rawIntColumnEnd <- rawIntColumnStart+nbSamples-1
spectralCountColumnsStart <-  rawIntColumnEnd +1
spectralCountColumnsEnd <- spectralCountColumnsStart+nbSamples-1			

if(opt$verbose){
	print("unique peptides after summation")
	print(length(rownames(dataSummed)))
}

### CREATE DATA FRAME OF UNIQUE PEPTIDES END

### FILTER PEPTIDE LIST (PROTEIN NAME, PEPTIDE COUNTS)

if(opt$verbose){
	print("FILTERING PEPTIDE LIST")
}

proteinNameFilterCond <- regexpr(selectedProteinName,dataSummed$proteins,ignore.case=TRUE) > -1

modifNameFilterCond <- rep(1,length(dataSummed$modif))
if(selectedModifName != "." ){
	modifNameFilterCond <- regexpr(selectedModifName,dataSummed$modif,ignore.case=TRUE) > -1
}

nonConProteinFilter  <- !isCon(dataSummed$proteins)
#peptidesForQuantFilterCond <- dataSummed$Peptides.used.for.quantitation >= peptidesForQuantCutoff

peptidesBeforeFilter <- nrow(dataSummed)
if(opt$verbose){
	print("REMOVED ENTRIES")
	print(dataSummed[!(proteinNameFilterCond & nonConProteinFilter & modifNameFilterCond ),]$proteins)
}
dataSummed <- dataSummed[proteinNameFilterCond & nonConProteinFilter & modifNameFilterCond,]

peptidesAfterFilter <- nrow(dataSummed)

if(opt$verbose){
	print("PROTEIN NAME FILTER")
	print(paste("    ", peptidesBeforeFilter - peptidesAfterFilter," peptides filtered out"))
}

### FILTER PEPTIDE LIST (PROTEIN NAME,PEPTIDE COUNTS) END

### FILTER SCORE, MASS DIFF 
## (NOT SURE IF THIS IS A GOOD IDEA)

decoyCond.nf <- isDecoy(dataSummed$proteins)
scores.nf <- dataSummed$maxScores
massDiff.nf <- dataSummed$minMassDiffs

#if(scoreMassDiffFilter){
#	
#	### thrs 70% score (replace with fdr?)
#	probs <- c(1:20)/20
#	scoreQuantiles <- quantile(dataSummed$maxScores,probs)
#	scoreThrs <- scoreQuantiles[6]
#	
#	### thrs 90% mass
#	massDiffQuantiles <- quantile(dataSummed$minMassDiffs,probs)
#	massDiffThrs <- massDiffQuantiles[length(massDiffQuantiles)-2]
#	
#	peptidesBeforeFilter <- nrow(dataSummed)
#	dataSummed <- dataSummed[!(dataSummed$maxScores < scoreThrs & dataSummed$minMassDiffs > massDiffThrs),]
#	peptidesAfterFilter <- nrow(dataSummed)
#	if(opt$verbose){
#		print("scoreMassDiffFilter")
#		print(paste("    ", peptidesBeforeFilter - peptidesAfterFilter," peptides filtered out"))
#	}
#}

peptidesBeforeFilter <- nrow(dataSummed)
	dataSummed <- dataSummed[!(dataSummed$minMassDiffs > precursorMassFilter),]
	peptidesAfterFilter <- nrow(dataSummed)
	if(opt$verbose){
		print("precursorMassFilter")
		print(paste("    ", peptidesBeforeFilter - peptidesAfterFilter," peptides filtered out"))
}

### FILTER SCORE, MASS DIFF END


### CALCULATE FDR

if(opt$verbose){
	print("---- FDR CALCULATION ----")
}

decoyCond <- isDecoy(dataSummed$proteins)
scores <- dataSummed$maxScores
decoyScores <- scores[decoyCond]
targetScores <- scores[!decoyCond]
decoyMassDiff <- dataSummed$minMassDiffs[decoyCond]
targetMassDiff <- dataSummed$minMassDiffs[!decoyCond]
qvals <- getQvals(scores,decoyCond)

### CALCULATE FDR END

### FILTER PEPTIDE LIST (FDR, DECOY)

if(opt$verbose){
	print("---- FILTER VALID PROTEIN IDS ----")
}

fdrFilterCond <- qvals < fdrCutoff

peptidesBeforeFilter <- nrow(dataSummed)
dataSummed <- dataSummed[!decoyCond & fdrFilterCond,]
peptidesAfterFilter <- nrow(dataSummed)
if(opt$verbose){
	print(paste("    ", peptidesBeforeFilter- peptidesAfterFilter," peptides filtered out"))
	print(paste("Tot valid peptides",peptidesAfterFilter))
}

### FILTER PEPTIDE LIST (FDR, DECOY) END

### EXTRACT RAW DATA

#protein.names <- dataSummed$proteins
rawDataColNames <- gsub(".1$","",paste("raw",names(dataSummed)[rawIntColumnStart:rawIntColumnEnd],sep=""))
rawData <- data.frame(dataSummed[rawIntColumnStart:rawIntColumnEnd])
names(rawData) <- rawDataColNames
#rownames(rawData) <- protein.names
rawData.log <- log(rawData)

### EXTRACT RAW DATA END

### NORMALIZE INTENSITIES

if(opt$verbose){
	print("NORMALIZE INTENSITIES")
}

nomalization.gainFactors  <- getNormalizationFactors(rawData)
normData <- normalizeIntensities(rawData,nomalization.gainFactors)

### NORMALIZE INTENSITIES END

### REPLACE ZEROS

if(opt$verbose){
	print("Replace zeros")
}

### @TODO TEMP
#### TEMP FIX, BETTER TO DERIVE FROM DISTRIBUTION !!!!! 
#minIntensity <- min(normData[normData > 0])
minIntensity.log <- log(minIntensity)
#normData[normData == 0] <- minIntensity
normData[normData < minIntensity] <- minIntensity
normData.log <- log(normData)

### REPLACE ZEROS END

### CREATE EXPRESSION DATASET
eset <- createExpressionDataset(normData.log,expDesign)
#eset <- createExpressionDataset(normData,expDesign)
e <- exprs(eset)

### CREATE EXPRESSION DATASET

if(opt$verbose){
	print("Outlier detection")
}

### OUTLIER DETECTION

conditionNas <- data.frame(row.names = rownames(rawData))
conditionSd <- data.frame(row.names = rownames(rawData))
conditionLargeSd <- data.frame(row.names = rownames(rawData))

i <- 1
for(condName in condition.names){
	
	cols <- which(eset$condition == condName)
	nas <- list()
	sd <- list()
	largeSd <-list()
	
	#### TODO 1REP
	if(length(cols) > 1){
		nas[[orgCondNames[i]]] <- apply(rawData[,cols], 1 ,FUN=min) < minIntensity 
		
		sd.temp <- apply(rawData[,cols], 1 ,FUN=sd)
		sd[[orgCondNames[i]]] <- sd.temp
		
		
		largeSd[[orgCondNames[i]]] <- sd.temp >= apply(rawData[,cols], 1 ,FUN=median)
		
	}else{
		nas[[orgCondNames[i]]] <- rawData[,cols]< minIntensity 
		sd.temp <- sd(rawData[,cols])
		sd[[orgCondNames[i]]] <- sd.temp
		largeSd[[orgCondNames[i]]] <- sd.temp >= median(rawData[,cols])
		
	}
	
	
	conditionNas <- data.frame(conditionNas,nas)
	conditionSd <- data.frame(conditionSd,sd)
	conditionLargeSd <- data.frame(conditionLargeSd,largeSd)
	i <- i+1
}	



###  OUTLIER DETECTION END


### RATIO ANALYSIS

#f <- factor(as.character(eset$condition))
#design <- model.matrix(~f)
#fit <- eBayes(lmFit(eset,design))
#
#### median
#conditionMedianIntensities <- calculateMedianIntensitiesPerCond(normData,expDesign)
#medianRatios <- calculateRatios(conditionMedianIntensities)
#
#### mean
#conditionMeanIntensities <- calculateMeanIntensitiesPerCond(normData,expDesign)
#meanRatios <- calculateRatios(conditionMeanIntensities)

### median
conditionMedianIntensities <- calculateMedianIntensitiesPerCond(normData,expDesign)
medianRatios <- calculateRatios(conditionMedianIntensities)
names(conditionMedianIntensities) <- orgCondNames
names(medianRatios) <- orgCondNames[2:length(orgCondNames)]

### mean
#conditionMeanIntensities <- calculateMeanIntensitiesPerCond(normData,expDesign)
#meanRatios <- calculateRatios(conditionMeanIntensities)

#### NO STATISTICAL ANALYSIS AND GRAPHICS IF ONLY ONE REPLICATE PER CONDITION
if(enoughReplicates) {
	
	f <- factor(as.character(eset$condition))
	design <- model.matrix(~f)
	fit <- eBayes(lmFit(eset,design))
	
}else{
	cat("WARNING: NOT ENOUGH REPLICATES PER CONDITION -> NO STAT ANALYSIS, TOP X OPTION DISABLED\n")
	### TURN OFF GRAPHICS
	
	topX <- 0
	volcanoPlots <- FALSE
	mvaPlots <- !regexpr("m",opt$selectedGraphics) > -1
	correlationPlots <- !regexpr("c",opt$selectedGraphics) > -1
	hClustPlot <- FALSE
	deFDRPlot <- FALSE
	
}

### RATIO ANALYSIS END

### CALCULATE SPECTRAL COUNTS PER CONDITION

spectralCounts <- spectralCountsPerCondition(dataSummed,expDesign,spectralCountColumnsStart, rownames(dataSummed))
names(spectralCounts) <- orgCondNames

### CALCULATE SPECTRAL COUNTS PER CONDITION END

###### OUTPUT #####################################################################
print("Preparing Graphics and csv Exports")
### GRAPHICS
if(opt$verbose){
	print("CREATING PLOTS")
}
pdf(pdfFile)

### COLORS
sample.names <- rownames(pData(eset))
RGBColors <- col2rgb(colors()[1:length(colors())])
HSVColors <- rgb2hsv( RGBColors[1,], RGBColors[2,], RGBColors[3,],
		maxColorValue=255)
HueOrder <- order( HSVColors[1,], HSVColors[2,], HSVColors[3,] )
myCol <- colors()[HueOrder]
tmpCol <- length(myCol)
sampleColors <- myCol[seq(2,tmpCol,round(tmpCol/length(sample.names)))]
#conditionColors <- myCol[seq(2,tmpCol,round(tmpCol/length(condition.names)))]
conditionColors <- 1:length(condition.names)
### COLORS END

### EXP DESIGN PLOT

if(dispExpDesign){
	
	xlim <- c(-1,6)
	ylim <- c(-2,nbSamples+2)
	
	plot(0,0,type="n", xlim=xlim, ylim=ylim, main="Experimental Design", axes=FALSE, xlab="", ylab="")
	
	condYPosStep <- (nbSamples+2)/(nbConditions+1)
	sampleNb <- 1
	
	for(condNb in 1:nbConditions){
		
		cond <- names(expDesign)[condNb]
		text(1,(condNb)*condYPosStep,orgCondNames[condNb], col=cond)
		
		for(i in 1:expDesign[[as.character(condNb)]]){
			text(4,sampleNb,sample.names[sampleNb], col=cond)
			sampleNb <- sampleNb + 1
		}
	}
}

### EXP DESIGN PLOT

### FDR PLOTS

if(fdrPlots){
	
#	par(mfrow=c(2,1))
#	plot(density(c(targetMassDiff,decoyMassDiff)))
#	plot(density(c(targetScores,decoyScores) ))
#	hist(c(targetMassDiff,decoyMassDiff),breaks=100)
#	hist(c(targetScores,decoyScores) ,breaks=100)
#	par(mfrow=c(1,1))
	
	### score mass diff. scatter plots
	plot(massDiff.nf, scores.nf 
			, col = decoyCond.nf +1
			, pch=20, ylab="score", xlab="mass diff. (ppm)")
	
#	if(scoreMassDiffFilter){
#			lines(c(massDiffThrs,max(massDiffQuantiles)),c(scoreThrs,scoreThrs),col="blue",lwd=2)
#			lines(c(massDiffThrs,massDiffThrs),c(min(scores.nf),scoreThrs),col="blue",lwd=2)
#	}
	lines(c(precursorMassFilter,precursorMassFilter),c(min(scores.nf),max(scores.nf)),col="blue",lwd=2)
	legend("topright",c("target","decoy"), pch=20, col= c(1,2))
	
	### score mass diff. scatter plots END
	
	par(mfrow=c(3,1))
	plotScoreDistrib(targetScores,decoyScores,nbBins=100,scoreName="Confidence score",ylab="Peptide Counts",title="")
	plotROC(qvals[!decoyCond],fdrMax=0.05,nbDataPoints=100,ylab="Valid Peptides")
	plot(sort(scores),rev(sort(qvals)), xlim=c(0,200), xlab="Confidence score", ylab="FDR", type="l")
	par(mfrow=c(1,1))

}

### FDR PLOTS END

### INTENSITY CORRELATION PLOTS 

if(correlationPlots){
	#source("/Users/erikahrne/dev/R/workspace/Lib/BasicPlots.R")
	pairs.annot(normData.log,diagText =nomalization.gainFactors,main="Corr")
	pairs.annot(log(conditionMedianIntensities))
	
}

### INTENSITY CORRELATION PLOTS END

### MVA PLOTS

if(mvaPlots){
	mva.pairs(e, cex=0.6, log=FALSE)
}

### MVA PLOTS END

### INTENSITY DENSITY PLOTS

if(intensityDensityPlots){
	
	### RAW DATA
	par(mfrow=c(2,1), cex=0.6)
	rawDataHistograms <- list()
	maxIntensity <- max(rawData.log[rawData.log > 0], na.rm=TRUE)[1] +1
	minIntensity <- min(rawData.log[rawData.log > 0], na.rm=TRUE)[1] -1
	
	maxHistCounts <- 0
	breaks <- seq(minIntensity,maxIntensity,0.4)
	for(i in 1:length(sample.names)){
		h <- hist(rawData.log[,i][rawData.log[,i] > 0],plot=FALSE, breaks=breaks)
		rawDataHistograms[[i]] <- h
		
		mhc <- max(h$counts)[1]
		if(mhc > maxHistCounts){
			maxHistCounts <- mhc 
		}
	}
	
	plot(rawDataHistograms[[1]]$mids,rawDataHistograms[[1]]$counts
		, type="n"
		, xlim=c(minIntensity,maxIntensity)
		,ylim=c(0,maxHistCounts)
		,xlab="raw intensity"
		,ylab="counts"
		, main= "Raw Intensity Distribution")
	for(i in 1:length(sample.names)){
		lines(rawDataHistograms[[i]]$mids,rawDataHistograms[[i]]$counts ,col=sampleColors[i])
	}
	
	### barplot of raw intensity sums per condition
	mp <- barplot(apply(rawData,2,FUN=function(t){return(sum(t,na.rm=TRUE) )}),
		, xaxt = "n"
		,  xlab = ""
		, las=2
		, col=sampleColors[1:(length(sample.names))]
		, main="raw intensity sum" 
	)
	
	axis(1, labels = FALSE)
	text(mp , par("usr")[3] -0.5, srt = 45, adj = 1,
				labels = sample.names, xpd = TRUE)
	
	par(mfrow=c(1,1), cex=1)
	
	### RAW DATA END
	
	### raw data
#	rawData.log.density <- density(rawData.log[rawData.log > 0])
#	xlim.raw= c(min(rawData.log.density$x)*0.95,max(rawData.log.density$x)*1.05)
#	
#	plot(density(rawData.log[,1]), type="n", main="Raw Data", xlab="Intensity",xlim=xlim.raw)
#	i <- 1
#	for(sample in sample.names){
#		lines(density(rawData.log[,i]),col=sampleColors[i])
#		
#		#print(paste(sample,mean(rawData.log[,i][rawData.log[,i]>0],na.rm=TRUE)))
#		
#		i <- i+1
#	}
#	legend("topright",sample.names,fill=sampleColors[1:(i-1)])
	
	### norm data
	normData.log.density <- density(normData.log)
	xlim.norm= c(min(normData.log.density$x)*0.95,max(normData.log.density$x)*1.05)
	
	density.sample1 <- density(normData.log[,1])
	ymax <- max(density.sample1$y)[1]*1.1
	
	naFraction <- sum(rawData[,1] == 0) / length(rawData[,1]) 
	textYPos <- max(density.sample1$y)* naFraction * 8
	plot(density.sample1, type="n", main="Normalized Data", xlab="Intensity", xlim=xlim.norm, ylim=c(0,ymax))
	points(minIntensity.log,textYPos, type="h" ,lty=3)
	text(minIntensity.log,textYPos,"Replaced\nNA's", pos=3,cex=0.9)
	i <- 1
	for(cond in sample.names){
		lines(density(normData.log[,i]),col=sampleColors[i])
		i <- i+1
	}
	legend("topright",sample.names,fill=sampleColors[1:(i-1)])
	
	
	if(enoughReplicates){
		### norm data per condition
		indexCond <- which(eset$condition == condition.names[1])
		plot(density.sample1, type="n", main="Normalized Data per Condition", xlab="Median Intensity", xlim=xlim.norm, ylim=c(0,ymax))
		points(minIntensity.log,textYPos, type="h" ,lty=3)
		text(minIntensity.log,textYPos,"Replaced\nNA's", pos=3,cex=0.9)
		#plot(density(rowMedians(e[,indexCond])), type="n", main= "Normalized Intensity per Condition", xlab="Median Intensity")
		i <- 1
		for(cond in condition.names){
			
			indexCond <- which(eset$condition == cond)
			lines(density(rowMedians(e[,indexCond])),col=conditionColors[i])
			i <- i+1
		}
		legend("topright",orgCondNames,fill=conditionColors[1:length(conditionColors)])
	}
}

### INTENSITY DENSITY PLOTS END

#### SPECTRAL COUNTS DENSITY PLOTS 
#if(spectralcountsDensityPlots){
#	### spectral count per coondition
#	
#	density.sample1 <- density(spectralCounts[,1])
#	ymax <- max(density.sample1$y)[1]*1.1
#	
#	plot(density.sample1, type="n", main= "Spectral Count", xlab="spectral counts per protein", ylim=c(0,ymax))
#	i <- 1
#	for(cond in condition.names){
#		lines(density(spectralCounts[,i]),col=conditionColors[i])
#		i <- i+1
#	}
#	legend("topright",orgCondNames,fill=conditionColors[1:length(conditionColors)])
#	
#}
#
#### SPECTRAL COUNTS DENSITY PLOTS END

### VOLCANO PLOTS

if(volcanoPlots){
	
	plotVolcanoes(fit, medianRatios,orgCondNames,conditionNas,conditionLargeSd,ratioCutOff=ratioCutOff
			,qvalVolcano=qvalVolcano, deFdrCutoff=deFdrCutoff,BCutoff=BCutoff)
	
}

#### VOLCANO PLOTS END

### HIERARCHICAL CLUSTERING

if(hClustPlot){
	
	if(FALSE){
		hClustHeatMap(eset, condNames=condition.names,  orgCondNames=orgCondNames)
	}
	
	
	### selection
	#selected.temp <-  (fit$lods[,2:ncol(fit)] < deFdrCutoff) & (abs(fit$coefficients[,2:ncol(fit)]) > log(ratioCutOff))
	selected.temp <-  (p.adjust(fit$p.value[,2:ncol(fit)],method="fdr") < deFdrCutoff) & (abs(fit$coefficients[,2:ncol(fit)]) > log(ratioCutOff))
	
	selected <- (apply(as.matrix(selected.temp), 1 ,FUN=sum) > 0)
	if(sum(selected) > 4){
		esetSel <- eset[selected,]
		hClustHeatMap(esetSel, main=paste("FDR Cut off: ", deFdrCutoff ,"& Ratio Cut off:", ratioCutOff ), condNames=condition.names, orgCondNames=orgCondNames)
	}else{
		print(paste("NO hClustPlot of selected features: Not enough entries meeting selection critera deFdrCutoff > ",deFdrCutoff," and ratioCutOff > ", ratioCutOff))
	}
}

### HIERARCHICAL CLUSTERING END

#### LODS DISTRUBUTION PER CONDITION
#if(lodsDistributionPlots){
#	plot(density(fit$lods[,2]),col=conditionColors[1], xlab="Log Odds",main="")
#	if(length(condition.names) >= 3){
#		for(i in 3:length(condition.names)){
#			
#			lines(density(fit$lods[,i]),col=conditionColors[i-1])
#		}
#	}
#	legend("topright",condition.names[2:length(condition.names)],fill=conditionColors[1:(length(condition.names)-1)])
#}
#
#### LODS DISTRUBUTION PER CONDITION END

### D.E. FDR PLOTS

if(deFDRPlot){
	
	plotDeFDR(fit,orgCondNames )
	
}

### D.E. FDR PLOTS END

##### TEMP P-VAL DISTRUBUITIONS
#
#for(i in 2:length(condition.names)){
#	hist(fit$p.value[,i], main=condition.names[i], breaks =100)
#}
#
##### TEMP P-VAL DISTRUBUITIONS END

cat(paste("Plots exported to",pdfFile),"\n")


##### TEMP @TODO
#p.val <- c()
#for(i in 1:nrow(normData)){
#	p.val <- c(p.val,t.test(normData[i,1:3],normData[i,4:6], alternative = c("two.sided"), var.equal = FALSE)$p.value)
#	#p.val <- c(p.val,t.test(log10(normData[i,1:3]),log10(normData[i,4:6]), alternative = c("two.sided"))$p.value)
#	
#} 
#
#par(mfrow=c(2,1))
#plot(log10(conditionMedianIntensities[,1]),log10(p.val))
#plot(log10(conditionMedianIntensities[,1]),log10(fit$p.value[,2]))
#
#par(mfrow=c(1,1))
#
##print(p.val)
#
#plot(p.val,fit$p.value[,2], main=round(cor(p.val,fit$p.value[,2]),2))
#points(p.val,fit$p.value[,2], col=2*(conditionLargeSd[,1]*1), pch=19)
#points(p.val,fit$p.value[,2], col=3*(conditionLargeSd[,2]*1), pch=19)
#points(p.val,fit$p.value[,2], col=4*(conditionLargeSd[,3]*1), pch=19)
#abline(lm(fit$p.value[,2] ~ p.val))
#
#plot(log2(medianRatios[,1]),-log10(p.val))
#plot(log2(medianRatios[,1]),-log10(p.adjust(p.val,method="fdr")))
#
#par(mfrow=c(2,1))
#
#hist(p.val, breaks = 100)
#plot(p.val,p.adjust(p.val,method="fdr"))
#
#
#for(i in 2:length(colnames(fit$p.value))){
#	hist(fit$p.value[,i], breaks = 100)
#	plot(fit$p.value[,i],p.adjust(fit$p.value[,i],method="fdr")) 
#}

##### TEMP END


graphics.off()
### GRAPHICS END

### CSV FILE

colnames(normData) <-  paste(colnames(normData) , "normInt", sep="")
names(conditionMedianIntensities) <- paste(names(conditionMedianIntensities), "medianInt", sep="")
names(medianRatios) <- paste(names(medianRatios),"vs",orgCondNames[1],"ratio", sep="")
#names(meanRatios) <- paste(names(meanRatios),"vs",condition.names[1],"meanRatio", sep="")

names(conditionNas) <- paste(names(conditionNas), "NA" , sep="")
names(conditionSd) <- paste(names(conditionSd), "Sd" , sep="")
names(conditionLargeSd) <- paste(names(conditionLargeSd), "LargeSd" , sep="")
names(spectralCounts) <- paste(names(spectralCounts), "SpecCount" , sep="")
### qval
#qvalExport <- data.frame(row.names = rownames(fit$p.value))
#for(i in 2:length(colnames(fit$p.value))){
#	qvalExport <- data.frame(qvalExport ,p.adjust(fit$p.value[,i],method="fdr")) 
#}
#colnames(qvalExport) <-  paste(condition.names[2:length(condition.names)] ,"vs",condition.names[1],"qval", sep="")
if(enoughReplicates){
	qvalExport <- data.frame(row.names = rownames(fit$p.value))
	for(i in 2:length(colnames(fit$p.value))){
		qvalExport <- data.frame(qvalExport ,p.adjust(fit$p.value[,i],method="fdr")) 
	}
	colnames(qvalExport) <-  paste(orgCondNames[2:length(condition.names)] ,"vs",orgCondNames[1],"qval", sep="")
}else{
	qvalExport <- rep("NA",length(dataSummed$peptides))
}

#dataSummed$peptides <- gsub(" ","",dataSummed$peptides)

results <- data.frame(dataSummed$peptides
		,dataSummed$proteins
		,dataSummed$descriptions
		,dataSummed$maxScores
		,dataSummed$modifs
		,dataSummed$minMassDiffs
		,dataSummed$missedCleavageSites
		,dataSummed$summedEntries
		,rawData
		,normData
		,conditionNas*1
		,conditionSd
		,conditionLargeSd*1
		,conditionMedianIntensities
		,medianRatios
		#,meanRatios
		,qvalExport
		,spectralCounts
		,row.names = row.names(dataSummed)) 

names(results) <- gsub("dataSummed\\.","",names(results))
names(results) <- gsub("s$","",names(results))
#write.matrix(results,file=csvFile, sep="$")
write.table(results,file=csvFile, sep="\t", row.names=FALSE)

### CSV FILE END

cat(paste("Results written to",csvFile),"\n")

### WRITE TOP X CSV

#TODO
# Give average intensity for top X for each condition.. indicate if fewer the X peptides were used.
# Note that a different set of peptides may be used for the same protein across different conditions.
#	protein, condYtopXAvgInt, condYNbPeptides, condYPeptideI, ... condYPeptideX, ....  	



if(topX > 0){
	
	if(opt$verbose){
		print(paste("WRITING TOP X FILE",topXCsvFile))
	}
	
	proteins <- dataSummed$proteins
	proteins.unique <- unique(proteins)
	
	header = c("protein"
				,paste(orgCondNames,"TopXMeanInt",sep="")
				,paste(orgCondNames,"TotNbPeptides",sep="") 
				,paste(orgCondNames,"TopXNbPeptides",sep="") 
				,paste(orgCondNames,"Modifs",sep="")
				)
		
	norm.tmp <- data.frame(normData)
					
	cat(c(paste(header,sep="\t"),"\n"),file=topXCsvFile,append=FALSE,sep="\t")
	for(protein in  proteins.unique){
		
		protein.indices =  which(protein == proteins)
				
		currentCondStartSample = 1
		meanInt = c()
		nbPeptides = c()
		topXNbPeptides = c()
		modifs = c()
		
		for(cm in condition.names){
			
			currentCondEndSample = currentCondStartSample + expDesign[[as.character(cm)]] -1
			
			selectedFrame.norm = norm.tmp[protein.indices,currentCondStartSample:currentCondEndSample]
			selectedFrame.raw = rawData[protein.indices,currentCondStartSample:currentCondEndSample]
			
			topXmedianInt.norm = apply(selectedFrame.norm,1, FUN=median)
			topXmedianInt.raw = apply(selectedFrame.raw,1, FUN=median)
			
			topXmedianInt.norm.order = order(topXmedianInt.norm, decreasing=TRUE)
		
			ordered.selected = c()
			if(length(topXmedianInt.norm.order) >= topX){
				ordered.selected = topXmedianInt.norm.order[1:topX]
			}else{
				ordered.selected = topXmedianInt.norm.order[1:length(topXmedianInt.norm.order)]
			}	
			
			meanInt = c(meanInt, mean(topXmedianInt.norm[ordered.selected], na.rm=TRUE))
			pepCount = sum(topXmedianInt.raw > 0 )
			nbPeptides = c(nbPeptides,pepCount)
			
			topXCount = topX
			if(pepCount < topX){
				topXCount = pepCount
			}
			
			topXNbPeptides = c(topXNbPeptides, topXCount)
			
			mofidTags = c()
			for(tmp in rownames(selectedFrame.norm)){
				
				tag = gsub("\\.\\.","",unlist(strsplit(tmp,"_"))[2])
				tag = gsub("^\\.$","",tag)
				modifTags = c(mofidTags,tag)
			}
			
			modifs = c(modifs,paste(modifTags,collapse="::"))
			
			currentCondStartSample = currentCondStartSample + expDesign[[as.character(cm)]] 
		}
		
		out <- data.frame(matrix(c(protein,meanInt,nbPeptides,topXNbPeptides, modifs), ncol=length(header)))
		write.table(out,file=topXCsvFile,append=TRUE,sep="\t", col.names=FALSE, row.names=FALSE)
		
	}
	
	cat(paste("Top X Results written to",topXCsvFile),"\n")
	
}

### WRITE TOP X CSV END

### save RObject

if(opt$saveRObject){
	
	rDataFile <- paste(outputDir,resultsFileLabel,".rData",sep="")
	save(results,eset,fit,condition.names, orgCondNames,file=rDataFile)
	cat(paste("RObject written to",rDataFile),"\n")
	
}

### save RObject end

###### OUTPUT END #################################################################

cat("DONE\n")
