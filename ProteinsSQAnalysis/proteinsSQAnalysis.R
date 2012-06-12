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
			
			make_option(c("-p", "--peptidesForQuantCutoff"), type="integer", default=1,
					help="Only include those proteins with have 
		at least x identified peptides [default %default]",
					metavar="number of peptides"),
			
			make_option(c("-c", "--controlCondition"), type="integer", default=1,
					help="Specifiy condition number to be used as control [default %default]",
					metavar="selected control condition"),
			
			
			make_option(c("-f", "--fdrCutoff"), type="double", default=0.01,
					help="Calculating intetnsity ratios for proteins with 
		a score corresponding to a fdr cut-off < specified value. [0-1]
		[default %default]",
					metavar="Protein level FDR cutoff"),
					
			make_option(c("-r", "--ratioCutoff"), type="double", default=1.5,
					help="Intensity ratio (fold change) cut-off used for
		graphics and results export. >1 [default %default]",
					metavar="Intensity ratio cutoff"),
			
			make_option(c("-q", "--deFdrCutoff"), type="double", default=0.01,
					help="fdr cut-off used for graphics. 
		High-lighting features with a qval < specified value. [0-1] [default %default]",
					metavar="Differential expression fdr cut-off"),
			make_option(c("-b", "--logOddsCutoff"), type="double", default=2.2,
					help="log odds cut-off used for graphics.
		High-lighting features with a log odds ratio > specified value. [default %default]",
					metavar="Intensity ratio logOdds cutoff"),
			make_option(c("-n", "--naReplace"), type="double", default=100,
					help="Intensity replacing missing values. [default %default]",
			),
			make_option(c("-x", "--proteinAccessionNormalization"), type="character", default=".",
					help="Normalize Intensities by selected protein(s)[default %default] (use all proteins).
		!!! Note that the specified proteins will be excluded from results output !!!",
					metavar="protein specific intensity normalization"),
	#		make_option(c("-m", "--meanFC"), action="store_true", default=FALSE,
	#				help="Display ratio of mean intensities per condition instead of median based intensity ratio [default %default]")
	#		,

			make_option(c("-s", "--saveRObject"), action="store_true", default=FALSE,
				help="Save R objects in 'label'.RData file [default %default]"),

			make_option(c("-z", "--selectedGraphics"), type="character", default="",
					help="Excluded Graphics: give letter for each plot to exclude ex: -z iv 
			(creates all plots but intensity density plots & volcano plot)
				experimental design (e)
				protein score distrib related plots (f)
				intensity density plots (i)
				volcano plots (v)
				sample and condition intensity correlation plots (c)
				sample MVA plots (m)
				hierarchical clustering plots (h)
				differential expression fdr plot (d)	
				qval based volcano plot instead of log odds based (l)
	 			[default (all plots) %default]")

)

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

normAC <- opt$proteinAccessionNormalization

peptidesForQuantCutoff <- opt$peptidesForQuantCutoff
if(is.na(opt$peptidesForQuantCutoff) | peptidesForQuantCutoff < 0 ){
	print(paste("ERROR. peptidesForQuantCutoff must be >= 0. You specified ", peptidesForQuantCutoff))
	q(status=-1)
}

###@TODO
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

minIntensity <- opt$naReplace
if(is.na(opt$naReplace) | minIntensity <= 0){
	print(paste("ERROR. deFdrCutoff must be > 0. You specified",minIntensity)) 
	q(status=-1)
}

#meanFC = opt$meanFC

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

#suppressPackageStartupMessages( library(MASS,quiet=TRUE))
suppressPackageStartupMessages(library(limma,quiet=TRUE))
#suppressPackageStartupMessages(library(Biobase,quiet=TRUE))
suppressPackageStartupMessages(library(affy,quiet=TRUE))
suppressPackageStartupMessages(library(gplots, quiet=TRUE ))
#source("/Users/erikahrne/dev/R/workspace/Lib/SearchResultsGraphics.R")
#source("/Users/erikahrne/dev/R/workspace/Lib/DecoyAnalysis.R")
#source(paste(getwd(),"/functions/ProgenesisPPFunctions.R",sep=""))

### get directory of script
initial.options <- commandArgs(trailingOnly = FALSE)
script.basename  <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
source(paste(script.basename,"/functions/ProgenesisPPFunctions.R",sep=""))

### LOAD SCRIPTS END

### PARSE CSV FILE

if(opt$verbose){
	print(paste("PARSING CSV FILE",progenesisFile))
}

	### READ EXPERIMENTAL DESIGN

line2header <- read.csv(file=progenesisFile,head=TRUE,sep=",",check.names=FALSE, skip=1,nrows=1, strip.white=TRUE)
line2header <- gsub(" ","",names(line2header))
line2header.names.unique <- unique(line2header)
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
for(cond in condition.names){
	expDesign[[cond]] <- 0
}

currentLab <- names(expDesign)[1] 
replicates <- 0

for(lab in line2header){
	
	if(lab %in% names(expDesign)){
		expDesign[[currentLab]] <- replicates + 1
		currentLab <- lab
		replicates <- 0
	}else{
		replicates <- replicates+1
	}
}

### get numbewr of samples
nbSamples <- 0
for(cond in condition.names){
	nbSamples <- nbSamples+ expDesign[[cond]]
}

enoughReplicates <- sum(unlist(expDesign)) >= 2*length(names(expDesign))

if(opt$verbose){
	print("Experimental design")
}	

	### READ EXPERIMENTAL DESIGN END


if(opt$verbose){
	print("condition.names")
}

### get number of conditions
nbConditions <- length(condition.names)

data <- read.csv(file=progenesisFile,head=TRUE,sep=",", skip=2)
col.names <- names(data)[regexpr("X",names(data)) == -1]

if(opt$verbose){
	print("col.names")
	print(col.names)
}

rawIntColumnStart <- 10+nbSamples
rawIntColumnEnd <- rawIntColumnStart+nbSamples-1
spectralCountColumnsStart <-  rawIntColumnEnd + 1

### PARSE CSV FILE END

### FILTER PROTEIN LIST (NAME, PEPTIDE COUNTS)

if(opt$verbose){
	print("FILTERING PROTEIN LIST")
}


proteinNameFilterCond <- regexpr(selectedProteinName,data$Accession,ignore.case=TRUE) > -1
nonConProteinFilter  <- !isCon(data$Accession)
peptidesForQuantFilterCond <- data$Peptides.used.for.quantitation >= peptidesForQuantCutoff

proteinsBeforeFilter <- nrow(data)
data <- data[proteinNameFilterCond & nonConProteinFilter & peptidesForQuantFilterCond,]
proteinsAfterFilter <- nrow(data)
if(opt$verbose){
	print(paste("    ", proteinsBeforeFilter- proteinsAfterFilter," proteins filtered out"))
}


### FILTER PROTEIN LIST (NAME,PEPTIDE COUNTS) END

### CALCULATE FDR

if(opt$verbose){
	print("---- FDR CALCULATION ----")
}

### ALL
decoyCond <- isDecoy(data$Accession)
scores <- data$Confidence.score
decoyScores <- scores[decoyCond]
targetScores <- scores[!decoyCond]
qvals <- getQvals(scores,decoyCond)


### CALCULATE FDR END

### FILTER PROTEIN LIST (FDR, DECOY)

if(opt$verbose){
	print("---- FILTER VALID PROTEIN IDS ----")
}

fdrFilterCond <- qvals < fdrCutoff

proteinsBeforeFilter <- nrow(data)
data <- data[!decoyCond & fdrFilterCond,]

proteinsAfterFilter <- nrow(data)
if(opt$verbose){
	print(paste("    ", proteinsBeforeFilter- proteinsAfterFilter," proteins filtered out"))
}

### FILTER PROTEIN LIST (FDR, DECOY) END

### EXTRACT RAW DATA

rawDataColNames <- gsub(".1$","",paste("raw",names(data)[rawIntColumnStart:rawIntColumnEnd],sep=""))
rawData <- data.frame(data[rawIntColumnStart:rawIntColumnEnd])
names(rawData) <- rawDataColNames
rownames(rawData) <- data$Accession

### EXTRACT RAW DATA END

### CALCULATE SPECTRAL COUNTS PER CONDITION

spectralCounts <- spectralCountsPerCondition(data,expDesign,spectralCountColumnsStart, data$Accession)
names(spectralCounts) <- orgCondNames

### CALCULATE SPECTRAL COUNTS PER CONDITION

### @TODO allow user to specify control condition
### reorder rawData, expDesign, orgCondNames, spectralCounts,  placing control condition in first position

if(selectedControlCond != 1){
	
	if(selectedControlCond <= nbConditions){
		rawData.tmp <-data.frame( row.names= rownames(rawData))
		rawData.control.tmp <- data.frame( row.names= rownames(rawData))
		
		expDesign.tmp <-list()
		
		orgCondNames.tmp <- c()
		
		spectralCounts.tmp <- data.frame( row.names= rownames(rawData))
		
		
		i <-1
		j <-2
		for(cond in names(expDesign)){
			
			condNbSamples <-  expDesign[[cond]]
			cond <- as.numeric(cond)
			
			if(cond == selectedControlCond){
				rawData.control.tmp <- rawData[,i:(i+condNbSamples-1)]
				expDesign.tmp[["1"]] <- condNbSamples
				
			}else{
				rawData.tmp <-data.frame(rawData.tmp, rawData[,i:(i+condNbSamples-1)] )
				expDesign.tmp[[as.character(j)]] <- condNbSamples
				orgCondNames.tmp <- c(orgCondNames.tmp, orgCondNames[cond])
				spectralCounts.tmp <-  data.frame(spectralCounts.tmp, spectralCounts[,cond])
				
				j <- j+1
			}  
			
			i <-i+condNbSamples
			
		}
		
		### put control cond first 
		rawData.tmp <- data.frame(rawData.control.tmp,rawData.tmp)
		orgCondNames.tmp <- c(orgCondNames[selectedControlCond],orgCondNames.tmp)
		
		spectralCounts.tmp <-  data.frame(spectralCounts[,selectedControlCond], spectralCounts.tmp)
		names(spectralCounts.tmp) <- orgCondNames.tmp
#	print("NEW")
#	print(names(rawData.tmp)) 
#	print(names(expDesign))
#	print(names(expDesign.tmp))
#	print(expDesign.tmp)
#	print(orgCondNames.tmp)
#	print(names(spectralCounts))
#	print(names(spectralCounts.tmp))
		
		### apply changes
		expDesign <- expDesign.tmp
		rawData <-  rawData.tmp
		orgCondNames <- orgCondNames.tmp
		spectralCounts <- spectralCounts.tmp 
	
	}else{
		
		print(paste("ERROR. Invalid controlCondition must be <=", nbConditions , ", You specified", selectedControlCond))
		q(status=-1)
		
	}
}





### NORMALIZE INTENSITIES

if(opt$verbose){
	print("NORMALIZE INTENSITIES")
}

normSelection <- regexpr(normAC,row.names(rawData),ignore.case=TRUE) > -1
		
if(sum(normSelection) == 0){
	print(paste("Error: Invalid protein selection for normalization",normAC))
	q(status=-1)
}		

nomalization.gainFactors  <- getNormalizationFactors(rawData[normSelection,])

if(opt$verbose){
	print("Normalization gain factors")
	print(nomalization.gainFactors)
}

### if normalization based on protein Ac filter out corresponding protein(s)
if(normAC != "."){
	rawData <- rawData[!normSelection,]
	data <- data[!normSelection,]
}

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
	
	#### @TODO 1REP
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

if(opt$verbose){
	print("Ratio Analysis")
}

### median
conditionMedianIntensities <- calculateMedianIntensitiesPerCond(normData,expDesign)
medianRatios <- calculateRatios(conditionMedianIntensities)
names(conditionMedianIntensities) <- orgCondNames
names(medianRatios) <- orgCondNames[2:nbConditions]

### mean
#conditionMeanIntensities <- calculateMeanIntensitiesPerCond(normData,expDesign)
#meanRatios <- calculateRatios(conditionMeanIntensities)

#### NO STATISTICAL ANALYSIS AND GRAPHICS IF ONLY ONE REPLICATE PER CONDITION
if(enoughReplicates) {

	f <- factor(as.character(eset$condition))
	
	design <- model.matrix(~f)
	
	fit <- eBayes(lmFit(eset,design))
	
	#print(topTable(fit, coef=3, adjust="fdr"))
	

}else{
	cat("WARNING: NOT ENOUGH REPLICATES PER CONDITION -> NO STAT ANALYSIS\n")
	### TURN OFF GRAPHICS
	
	volcanoPlots <- FALSE
	mvaPlots <- !regexpr("m",opt$selectedGraphics) > -1
	correlationPlots <- !regexpr("c",opt$selectedGraphics) > -1
	hClustPlot <- FALSE
	deFDRPlot <- FALSE
	
}

### RATIO ANALYSIS END



###### OUTPUT #####################################################################

### GRAPHICS
if(opt$verbose){
	cat("CREATING PLOTS\n")
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
	
	### ALL
	par(mfrow=c(3,1))
	plotScoreDistrib(targetScores,decoyScores,nbBins=100,scoreName="Confidence score",ylab="Protein Counts",title="")
	#plotScoreDistrib(targetScores,decoyScores,nbBins=100,scoreName="Confidence score",ylab="Protein Counts",title="", xlim=c(0,200))
	plotROC(qvals[!decoyCond],fdrMax=0.05,nbDataPoints=100,ylab="Valid Proteins")
	#plot(sort(scores),rev(sort(qvals)), xlim=c(0,max(scores)), xlab="Confidence score", ylab="FDR", type="l")
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
	
	rawData.log <- log(rawData)
	
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
	
	rawData.log <- log(rawData)
	### raw data
#	rawData.log.density <- density(rawData.log[rawData.log > 0])
#	xlim.raw= c(min(rawData.log.density$x)*0.95,max(rawData.log.density$x)*1.05)
#		
#	plot(density(rawData.log[,1]), type="n", main="Raw Data", xlab="Intensity",xlim=xlim.raw)
#	i <- 1
#	for(cond in sample.names){
#		lines(density(rawData.log[,i]),col=sampleColors[i])
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
			
	plot(density.sample1, type="n", main="Normalized Data", xlab="Intensity", xlim=xlim.norm , ylim=c(0,ymax))
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
		#indexCond <- which(eset$condition == condition.names[1])
		plot(density.sample1, type="n", main="Normalized Data per Condition", xlab="Median Intensity",xlim=xlim.norm , ylim=c(0,ymax))
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

### SPECTRAL COUNTS DENSITY PLOTS 
#if(spectralcountsDensityPlots){
#	### spectral count per coondition
#	plot(density(spectralCounts[,1]), type="n", main= "Spectral Count", xlab="spectral counts per protein")
#	i <- 1
#	for(cond in condition.names){
#		lines(density(spectralCounts[,i]),col=conditionColors[i])
#		i <- i+1
#	}
#	legend("topright",orgCondNames,fill=conditionColors[1:length(conditionColors)])
#	
#}

### VOLCANO PLOTS

if(volcanoPlots){
	
	plotVolcanoes(fit, medianRatios,orgCondNames,conditionNas,conditionLargeSd,ratioCutOff=ratioCutOff
		,qvalVolcano=qvalVolcano, deFdrCutoff=deFdrCutoff,BCutoff=BCutoff)
	
}

### VOLCANO PLOTS

### HIERARCHICAL CLUSTERING

if(hClustPlot){
	
	#hClustHeatMap(eset, condNames=orgCondNames)
	hClustHeatMap(eset, condNames=condition.names, orgCondNames)

	### selection
	#selected.temp <-  (fit$lods[,2:ncol(fit)] < deFdrCutoff) & (abs(fit$coefficients[,2:ncol(fit)]) > log(ratioCutOff))
	selected.temp <-  (p.adjust(fit$p.value[,2:ncol(fit)],method="fdr") < deFdrCutoff) & (abs(fit$coefficients[,2:ncol(fit)]) > log(ratioCutOff))
	
	selected <- (apply(as.matrix(selected.temp), 1 ,FUN=sum) > 0)
	if(sum(selected) > 4){
		esetSel <- eset[selected,]
		hClustHeatMap(esetSel, main=paste("FDR Cut off: ", deFdrCutoff ,"& Ratio Cut off:", ratioCutOff ), condNames=condition.names, orgCondNames)
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

cat(paste("Plots exported to",pdfFile),"\n")
graphics.off()

### GRAPHICS END

### CSV FILE

colnames(normData) <-  paste(colnames(normData) , "normInt", sep="")
names(conditionMedianIntensities) <- paste(orgCondNames, "medianInt", sep="")
#names(medianRatios) <- paste(condition.names,"vs",condition.names[1],"ratio", sep="")
names(medianRatios) <- paste(orgCondNames[2:nbConditions],"vs",orgCondNames[1],"ratio", sep="")

#names(meanRatios) <- paste(names(meanRatios),"vs",condition.names[1],"meanRatio", sep="")

names(conditionNas) <- paste(names(conditionNas), "NA" , sep="")
names(conditionSd) <- paste(names(conditionSd), "Sd" , sep="")
names(conditionLargeSd) <- paste(names(conditionLargeSd), "LargeSd" , sep="")
names(spectralCounts) <- paste(names(spectralCounts), "SpecCount" , sep="")

### qval
if(enoughReplicates){
	qvalExport <- data.frame(row.names = rownames(fit$p.value))
	for(i in 2:length(colnames(fit$p.value))){
		qvalExport <- data.frame(qvalExport ,p.adjust(fit$p.value[,i],method="fdr")) 
	}
	colnames(qvalExport) <-  paste(orgCondNames[2:nbConditions] ,"vs",orgCondNames[1],"qval", sep="")
}else{
	qvalExport <- rep("NA",length(data$Accession))
}



results <- data.frame(data$Accession
		,data$Peptide.count
		,data$Peptides.used.for.quantitation
		,data$Confidence.score
		,data$Description
		,rawData
		,normData
		,conditionNas*1
		,conditionSd
		,conditionLargeSd*1
		,conditionMedianIntensities
		,medianRatios
#		,meanRatios
		,qvalExport
		,spectralCounts
		,row.names = row.names(data)) 

names(results) <- gsub("data\\.","",names(results))
#write.matrix(results,file=csvFile, sep="$")
write.table(results,file=csvFile, sep="\t", row.names=FALSE)

### CSV FILE END

cat(paste("Results written to",csvFile),"\n")

### save RObject

if(opt$saveRObject){
	
	rDataFile <- paste(outputDir,resultsFileLabel,".rData",sep="")
	save(medianRatios,eset,fit,condition.names, orgCondNames,file=rDataFile)
	cat(paste("RObject written to",rDataFile),"\n")
	
}

### save RObject END 

cat("DONE\n")

###### OUTPUT END #################################################################

##### TEMP P-VAL DISTRUBUITIONS
#
#for(i in 2:length(condition.names)){
#	hist(fit$p.value[,i], main=condition.names[i], breaks =100)
#}
#
##### TEMP P-VAL DISTRUBUITIONS END

##### K-MEANS CLUSTERING
#
#par(mfrow=c(1,1))
#
#require(graphics)
#
## a 2-dimensional example
#
#x <- matrix(c(rnorm(99, mean = -10, sd = 0.05),rnorm(66, mean = 0, sd = 0.05),rnorm(99, mean = 10, sd = 0.05)),ncol=3)
#
#head(x)
#
##x <- rbind(matrix(rnorm(99, sd = 0.3), ncol = 3),
##		matrix(rnorm(99, mean = 1, sd = 0.3), ncol = 3),
##		matrix(rnorm(99, mean = 2, sd = 0.3), ncol = 3))
#
#colnames(x) <- c("x", "y")
#(cl <- kmeans(x, 3))
#plot(x, col = cl$cluster)
#points(cl$centers, col = 1:3, pch = 8, cex=2)
#
#print(cl$centers)
#
#kmeans(x,1)$withinss # if you are interested in that
#
### random starts do help here with too many clusters
#(cl <- kmeans(x, 5, nstart = 25))
#plot(x, col = cl$cluster)
#points(cl$centers, col = 1:5, pch = 8)
#
#
##### K-MEANS CLUSTERING END
