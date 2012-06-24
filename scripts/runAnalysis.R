# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

# 
# Author: erikahrne
###############################################################################
#### DEPENDANCIES
#limma
#affy
#gplots
#optparse

#####
#1) READ CMD LINE PARAMETERS
#2) PARSE CSV FILE		
#	a) PARSE prog DATA
#		- PARSE CSV	
#		- APPLY DATA FILTERS
#		- extract exp. design
#3) CREATE NORM DATA
#       - normalize to first sample
#		- replace missing values
#4) ANALYSIS OFDIFFERENTIAL EXPRESSION
#		- a) calculate median ratios
#		- b) calculate ratio sign. (qvalues)
#		- c) calculate cvs per cond
#5 ) CREATE GRAPHICS EXPORT
#6 ) CREATE .CSV EXPORT 

###################################### INIT

### SET IN MAIN SCRIPTS
# BASEDIR
# ISPEPTIDEANALYSIS

### SOURCE GENERAL FUNCTIONS
#### LOAD IDENTIFICATION RELATED FUNCTIONS (used by CsvParser)
source(paste(BASEDIR,"functions/IdentificationAnalysis.R", sep=""))
### LOAD PARSE FUNCTIONS
source(paste(BASEDIR,"functions/CsvParser.R", sep=""))
### LOAD EXPRESSION ANALYSIS FUNCTIONS
source(paste(BASEDIR,"functions/ExpressionAnalysis.R", sep=""))
### LOAD PLOTTING FUNCTIONS
source(paste(BASEDIR,"functions/Graphics.R", sep=""))


### GET VERSION
source(paste(BASEDIR,"config/VERSION.R", sep=""))
### END GET VERSION

if(ISPEPTIDEANALYSIS){
	
	######################################  CMD OPTIONS
	source(paste(BASEDIR,"PeptidesSQAnalysis/functions/CommandLineOptionsPeptides.R", sep=""))
	userOptions <- getUserOptions(version=VERSION)
	######################################  CMD OPTIONS END
	
	######################################  PARSE CSV FILE
	source(paste(BASEDIR,"PeptidesSQAnalysis/subscript/parsePeptides.R", sep=""))
	###################################### PARSE CSV FILE END
}else{
	
	######################################  CMD OPTIONS
	source(paste(BASEDIR,"ProteinsSQAnalysis/functions/CommandLineOptionsProteins.R", sep=""))
	userOptions <- getUserOptions(version=VERSION)
	######################################  CMD OPTIONS END
	
	######################################  PARSE CSV FILE
	source(paste(BASEDIR,"ProteinsSQAnalysis/subscript/parseProteins.R", sep=""))
	###################################### PARSE CSV FILE END
	
}

### CREATED DATA STRUCTURES
#PROGOUTDATA
#EXPDESIGN
#NBSAMPLES
#ISENOUGHREPLICATES
#NBCONDITIONS
#UNGROUPEDRAWINTDATA
#UNGROUPEDSPECCOUNTDATA

###################################### INIT END

###################################### EXPRESSION ANALYSIS
cat("EXPRESSION ANALYIS\n")
#### CREATE NORMALIZED AUC DATA

### LOAD EXPRESSION ANALYSIS FUNCTIONS
source(paste(BASEDIR,"functions/ExpressionAnalysis.R", sep=""))
### LOAD EXPRESSION ANALYSIS FUNCTIONS END

### normalize intensites
# create data structures UNGROUPEDNORMINTDATA, GROUPEDNORMINTDATA
if(userOptions$verbose){
	print("NORMALIZE INTENSITIES")
}
UNGROUPEDNORMINTDATA <- createUngroupedNormData(UNGROUPEDRAWINTDATA ,selectedNormAC=userOptions$normAC, userSpecifiedMinInt=userOptions$minIntensity,verbose=userOptions$verbose)

### group norm data per condition
GROUPEDNORMINTDATA <- getCondGroupedSampleData(UNGROUPEDNORMINTDATA,expDesign=EXPDESIGN)

#### CREATE NORMALIZED AUC DATA END

### STAT. ANALYSIS, CALCULATE RATIOS AND QVALUES OF DE FEATURES, CVS MARK CONDITIONS CONTAINING MISSING VALUES VARIATION

CONTROLCONDITION <- names(EXPDESIGN)[userOptions$selectedControlCond]
NONCONTROLCONDITIONS <-getNonControlConditionNbs(GROUPEDNORMINTDATA, controlCondition=CONTROLCONDITION)

### create data frame of ratios of medians
MEDIANNORMINTDATA <- getMedianIdxPerCond(GROUPEDNORMINTDATA)
### create data frame of median ratios per condition
RATIOSMEDIANNORMINTPERCOND <-  getMedianRatiosPerCondition(MEDIANNORMINTDATA
		,controlCondition=CONTROLCONDITION
		,nonControlConditions=NONCONTROLCONDITIONS)

### create data frame of DE qvalues
if(ISENOUGHREPLICATES) { # only if each condition has at least 2 replicates
	
	### LOAD EXT LIBRARIES
	suppressPackageStartupMessages(library(limma,quiet=TRUE))
	suppressPackageStartupMessages(library(affy,quiet=TRUE))
	suppressPackageStartupMessages(library(gplots, quiet=TRUE ))
	### LOAD EXT LIBRARIES END
	
#	QVALUESPERCOND <- getDEQvaluesPerCondition(GROUPEDNORMINTDATA
#			,controlCondition=CONTROLCONDITION
#			,nonControlConditions=NONCONTROLCONDITIONS
#			,expDesign=EXPDESIGN  )
	
	
	QVALUESPERCOND2 <- getPairWiseDEQvaluesPerCondition(GROUPEDNORMINTDATA
			, controlCondition=CONTROLCONDITION
			, nonControlConditions=NONCONTROLCONDITIONS
			, expDesign=EXPDESIGN)
	
	### calculate cv per condition
	CVSPERCOND <- getCVPerCondition(GROUPEDNORMINTDATA)
	
	### mark conditions with missing values (na)
#	NAPERCOND <- getNAPerCondition(getCondGroupedSampleData(UNGROUPEDRAWINTDATA,expDesign=EXPDESIGN)
#					,minIntensity=userOptions$minIntensity)
	
}else{
	cat("WARNING: NOT ENOUGH REPLICATES PER CONDITION -> NO QVALUES CALCULATED\n")
	### TURN OFF GRAPHICS
	userOptions$isVolcanoPlots <- FALSE
	userOptions$isCorrelationPlots <- FALSE
	userOptions$isHClustPlot <- FALSE
	userOptions$isDeFdrPlot <- FALSE
}

###  STAT. ANALYSIS, CALCULATE RATIOS AND QVALUES OF DE FEATURES, CVS, MARK CONDITIONS CONTAINING MISSING VALUES END

###################################### EXPRESSION ANALYSIS END


###### OUTPUT #####################################################################

######################################## GRAPHICS
cat("CREATING OUTPUT FILES\n")
pdf(userOptions$pdfFilePath)

if(userOptions$verbose){
	cat("CREATING PLOTS\n")
}

### COLORS
CONDITIONCOLORS <-  data.frame(COLORS[(1:NBCONDITIONS)], row.names=names(EXPDESIGN))
names(CONDITIONCOLORS) <- c("cond")

#SAMPLECOLORS <-  data.frame(1:length(names(UNGROUPEDRAWINTDATA)), row.names=names(UNGROUPEDRAWINTDATA))
### COLORS END

### EXP DESIGN PLOT
if(userOptions$isDispExpDesign){
	plotExpDesign(EXPDESIGN,sampleNames=names(UNGROUPEDRAWINTDATA), controlCondition=CONTROLCONDITION,condColors=CONDITIONCOLORS, version=VERSION )
}
### EXP DESIGN PLOT END

### FDR PLOTS
if(userOptions$isFdrPlots){
	
	if(ISPEPTIDEANALYSIS){
		### plot precursor mass diff vs. score
		plotMassDiffVsScore(FILTEROBJ$massDiffPlotObj$minMassDiffs,FILTEROBJ$massDiffPlotObj$scores,FILTEROBJ$massDiffPlotObj$isDecoy, precursorMassFilter=userOptions$precursorMassFilter)
	}

	### id related plots
	par(mfrow=c(3,1))
	plotScoreDistrib(FILTEROBJ$scores[!FILTEROBJ$decoyCond],FILTEROBJ$scores[FILTEROBJ$decoyCond],nbBins=100,scoreName="Confidence score",ylab="Protein Counts",title="")
	plotROC(FILTEROBJ$qvals[!FILTEROBJ$decoyCond],fdrMax=0.05,nbDataPoints=100,ylab="Valid Proteins")
	plotIdScoreVsFDR(FILTEROBJ$scores,FILTEROBJ$qvals, xlim=c(0,max(FILTEROBJ$scores)), xlab="Confidence score", ylab="FDR", type="l")
	par(mfrow=c(1,1))
}
### FDR PLOTS END

### INTENSITY CORRELATION PLOTS 
if(userOptions$isCorrelationPlots){
	
	### plots all against all per condition
	for(cond in names(EXPDESIGN)){
		pairs.annot(log(GROUPEDNORMINTDATA[[cond]]), main=cond)
	}
	
	### plot all against all condition median intensities
	pairs.annot(log(MEDIANNORMINTDATA))
	
	### boxplot of CVS per condition
	boxplot(data.frame(CVSPERCOND[names(EXPDESIGN)])*100,ylab= "C.V. (%)" ,col=as.character(CONDITIONCOLORS[names(EXPDESIGN),]),cex.lab=1, las=2 )
}

### INTENSITY CORRELATION PLOTS END

### INTENSITY DISTRIBUTION PLOTS
if(userOptions$isIntensityDensityPlots){
	
	### two plots per slide
	par(mfrow=c(2,1))
	plotIntensityDistributions(log10(UNGROUPEDRAWINTDATA)
			,isLegend=F
			, xlab="RAW INTENSITY"
			, ylab="FREQUENCY"
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGN)
			, lwd=1.5)
	
	sampleIntSumBarplot(UNGROUPEDRAWINTDATA
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGN)
			, main="RAW INT SUM PER SAMPLE")
	par(mfrow=c(1,1))
	
	plotIntensityDistributions(log10(UNGROUPEDNORMINTDATA)
			, isLegend=T
			, naReplacedInt=log10(userOptions$minIntensity)
			, xlab="NORM INTENSITY"
			, ylab="FREQUENCY"
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGN)
			, lwd=2)
	
	plotIntensityDistributions(log10(MEDIANNORMINTDATA)
			, isLegend=T
			, naReplacedInt=log10(userOptions$minIntensity)
			, xlab="PROTEIN COND. MEDIAN NORM INTENSITY"
			, ylab="FREQUENCY"
			, colors=as.character(CONDITIONCOLORS[,1])
			, lwd=2)
	
}

### INTENSITY DISTRIBUTION PLOTS END

### VOLCANO PLOTS
if(userOptions$isVolcanoPlots){
	plotAllVolcanoes(RATIOSMEDIANNORMINTPERCOND
			, QVALUESPERCOND2
			, CVSPERCOND
			, controlCondition=CONTROLCONDITION
			, nonControlConditions=NONCONTROLCONDITIONS
			, qvalueCutOff=userOptions$deFdrCutoff
			, ratioCutOff=userOptions$ratioCutOff)
}
### VOLCANO PLOTS END


### HIERARCHICAL CLUSTERING
if(userOptions$isHClustPlot){
	hClustHeatMap(UNGROUPEDNORMINTDATA
			, expDesign=EXPDESIGN
			, conditionColors=CONDITIONCOLORS )
}
### HIERARCHICAL CLUSTERING END

### D.E. FDR PLOTS
if(userOptions$isDeFdrPlot){
	
	### two plots per slide	
	par(mfrow=c(1,2))
	plotNbValidDeFeaturesPerFDR(QVALUESPERCOND2
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=TRUE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalCutOffs=seq(0,0.3, length.out=10)
			,conditionColors= CONDITIONCOLORS
			, main="UP-REG. PROTEINS"
    		,ylab="# PROTEINS"
	)
	
	plotNbValidDeFeaturesPerFDR(QVALUESPERCOND2
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=FALSE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalCutOffs=seq(0,0.3, length.out=10)
			,conditionColors= CONDITIONCOLORS	
			, main="DOWN-REG. PROTEINS"
			,ylab="# PROTEINS"
			,isLegend=F
	)
	par(mfrow=c(1,1))
}

### D.E. FDR PLOTS END
print(paste("CREATED FILE", userOptions$pdfFilePath))
graphics.off()
###################################### GRAPHICS END


###################################### TSV FILE EXPORT

### raw Intensities
rawDataExport <- UNGROUPEDRAWINTDATA
names(rawDataExport) <- paste("rawInt_",names(rawDataExport),sep="")

### norm Intensities
normDataExport <- UNGROUPEDNORMINTDATA
#names(normDataExport) <- paste("normInt_",names(normDataExport),sep="")

### median intensites 
medianNormIntExport <- data.frame(MEDIANNORMINTDATA)
names(medianNormIntExport) <- paste("medianNormInt_",names(MEDIANNORMINTDATA),sep="")

### spec count sums
specCountSumExport <- data.frame(getSummedIdxPerCond(getCondGroupedSampleData(UNGROUPEDSPECCOUNTDATA,expDesign=EXPDESIGN)))
names(specCountSumExport) <- paste("specCountSum_",names(EXPDESIGN),sep="")

### median ratios
medianRatiosExport <- data.frame(RATIOSMEDIANNORMINTPERCOND)
names(medianRatiosExport) <- paste("medianRatio_",names(RATIOSMEDIANNORMINTPERCOND),"_vs_",CONTROLCONDITION,sep="")

if(ISPEPTIDEANALYSIS){
	
	export <- data.frame(PROGOUTDATA$peptides
			,PROGOUTDATA$modifs
			,PROGOUTDATA$Accession
			,PROGOUTDATA$Description
			,PROGOUTDATA$Confidence.score
			,rawDataExport
			,normDataExport
			,medianNormIntExport
			,specCountSumExport
			,medianRatiosExport
	)
}else{
	
	export <- data.frame(PROGOUTDATA$Accession
			,PROGOUTDATA$Description
			,PROGOUTDATA$Peptides.used.for.quantitation
			,PROGOUTDATA$Confidence.score
			,rawDataExport
			,normDataExport
			,medianNormIntExport
			,specCountSumExport
			,medianRatiosExport
	)
}



### qvalues and coef. variance not calculated if only one sample per condition 
if(ISENOUGHREPLICATES){
	
	### d.e qvalues
	qvaluesExport <- data.frame(QVALUESPERCOND2)
	names(qvaluesExport) <- paste("qvalue_",names(QVALUESPERCOND2),"_vs_",CONTROLCONDITION,sep="")
	
	### coefficient of variance (%)
	cvsExport <- data.frame(CVSPERCOND*100)
	names(cvsExport) <- paste("cv_",names(CVSPERCOND),sep="")
	
	### add to export
	export <- data.frame(export,qvaluesExport,cvsExport)
}

### clean up column names
names(export) <- gsub("PROGOUTDATA.","",names(export))

write.table(export,file=userOptions$csvFilePath, sep="\t", row.names=FALSE)

print(paste("CREATED FILE", userOptions$csvFilePath))

###################################### TSV FILE EXPORT END


###################################### SAVE SERIALIZED R-DATA OBJECT

if(userOptions$isSaveRObject){
	
	rDataFile <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,".rData",sep="")
	
	save(GROUPEDNORMINTDATA
			,MEDIANNORMINTDATA
			,RATIOSMEDIANNORMINTPERCOND
			,QVALUESPERCOND2
			,CVSPERCOND		
			,CONTROLCONDITION
			,EXPDESIGN		
			,file=rDataFile)
	print(paste("CREATED FILE", rDataFile))
}

###################################### SAVE SERIALIZED R-DATA OBJECT END 

###### OUTPUT END #####################################################################

cat("DONE\n")










