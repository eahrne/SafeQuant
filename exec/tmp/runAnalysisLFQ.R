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
#RSVGTipsDevice

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
#6 ) CREATE .TSV EXPORT 

###################################### INIT
	
	### SET IN MAIN SCRIPTS
	# BASEDIR
	# ISPEPTIDEANALYSIS
	
	### SOURCE GENERAL FUNCTIONS
	#### LOAD IDENTIFICATION RELATED FUNCTIONS (used by CsvParser)
	source(paste(BASEDIR,"R/IdentificationAnalysis.R", sep=""))
	### LOAD PARSE FUNCTIONS
	source(paste(BASEDIR,"R/LFQCsvParser.R", sep=""))
	### LOAD EXPRESSION ANALYSIS FUNCTIONS
	source(paste(BASEDIR,"R/ExpressionAnalysis.R", sep=""))
	### LOAD PLOTTING FUNCTIONS
	source(paste(BASEDIR,"R/Graphics.R", sep=""))
	### LOAD SVG PLOTTING FUNCTIONS
	source(paste(BASEDIR,"R/SVGPlots.R", sep=""))
	
	### GET VERSION
	source(paste(BASEDIR,"config/VERSION.R", sep=""))
	### END GET VERSION
	
	source(paste(BASEDIR,"R/userOptions.R", sep=""))
	userOptions <- getUserOptions(version=VERSION)
	
	if(ISPEPTIDEANALYSIS){
		
		######################################  CMD OPTIONS END
		
		######################################  PARSE CSV FILE
		source(paste(LFQDIR,"PeptideAnalysis/subscript/parsePeptides.R", sep=""))
		###################################### PARSE CSV FILE END
	}else{
		
		######################################  CMD OPTIONS END
		
		######################################  PARSE CSV FILE
		source(paste(LFQDIR,"ProteinAnalysis/subscript/parseProteins.R", sep=""))
		###################################### PARSE CSV FILE END
		
	}
	
	### CREATED DATA STRUCTURES
	#PROGOUTDATA
	#EXPDESIGNOBJ
	#UNGROUPEDRAWINTDATA
	#UNGROUPEDSPECCOUNTDATA
	
	### DISCARD SAMPLES
	
	#test discard samples
	#userOptions$discardedSamples <-c(3)
	#discardedSamples <-c(6:19)
	
	if(length(userOptions$discardedSamples) > 0){
		
		if(userOptions$verbose){
			cat("DISCARDING SAMPLE(S) ",userOptions$discardedSamples ,"\n" )
		}
		
		EXPDESIGNOBJ <- discardSamplesFromExpDesignObj(EXPDESIGNOBJ,userOptions$discardedSamples )
		
		### remove samples from UNGROUPEDRAWINTDATA
		UNGROUPEDRAWINTDATA <- UNGROUPEDRAWINTDATA[ , !((1:ncol(UNGROUPEDRAWINTDATA))  %in% userOptions$discardedSamples  ) ]
		
		### remove samples from UNGROUPEDSPECCOUNTDATA
		UNGROUPEDSPECCOUNTDATA <- UNGROUPEDSPECCOUNTDATA[,   !((1:ncol(UNGROUPEDSPECCOUNTDATA))  %in% userOptions$discardedSamples  )  ]
		
	}
	
	### test user specified control condition
	if(userOptions$selectedControlCond >  EXPDESIGNOBJ$nbConditions){
		print(paste("ERROR: Invalid control condition specified ->",userOptions$selectedControlCond  ))
		quit(status=-1)
	}
	
	### DISCARD SAMPLES END
	
	### CHECK IF MORE THAN ONE CONDITION WITH MORE THAN ONE REP (IF NOT NO STAT ANALYSIS AND ACCOMPANYING POTS)
	ISDEEXP <- (min(EXPDESIGNOBJ$expDesign) > 1) & (ncol(EXPDESIGNOBJ$expDesign) > 1) ### @TODO get rid off this, add test to each function
	
	### SUPRESS WARNINGS
	if(!userOptions$verbose){
		options(warn=-1)
	}


###################################### INIT END

###################################### EXPRESSION ANALYSIS
	
	if(userOptions$verbose){
		print("EXP. DESIGN")
		print(EXPDESIGNOBJ)
	}
	
	cat("EXPRESSION ANALYIS\n")
	#### CREATE NORMALIZED AUC DATA
	
	### LOAD EXPRESSION ANALYSIS FUNCTIONS
	source(paste(BASEDIR,"R/ExpressionAnalysis.R", sep=""))
	### LOAD EXPRESSION ANALYSIS FUNCTIONS END
	
	### normalize intensites
	# create data structures UNGROUPEDNORMINTDATA, GROUPEDNORMINTDATA
	if(userOptions$verbose){
		print("NORMALIZE INTENSITIES")
	}
	
	### baselinInt -> int at 4 sd below mean
	BASELINEINTENSITY <- getBaselineIntensity(UNGROUPEDRAWINTDATA[,1])
	
	UNGROUPEDNORMINTDATA <- .createUngroupedNormData(UNGROUPEDRAWINTDATA ,selectedNormAC=userOptions$normAC, baselineIntensity=BASELINEINTENSITY,verbose=userOptions$verbose)
	
	### group norm data per condition
	GROUPEDNORMINTDATA <- getCondGroupedSampleData(UNGROUPEDNORMINTDATA,expDesign=EXPDESIGNOBJ$expDesign )
	
	#### CREATE NORMALIZED AUC DATA END
	
	### STAT. ANALYSIS, CALCULATE RATIOS AND QVALUES OF DE FEATURES, CVS MARK CONDITIONS CONTAINING MISSING VALUES VARIATION
	
	CONTROLCONDITION <- names(EXPDESIGNOBJ$expDesign )[userOptions$selectedControlCond]
	NONCONTROLCONDITIONS <-getNonControlConditionNames(GROUPEDNORMINTDATA, controlCondition=CONTROLCONDITION)
	
	### create data frame of ratios of medians
	MEDIANNORMINTDATA <- getMedianIdxPerCond(GROUPEDNORMINTDATA)
	### create data frame of median ratios per condition
	### Breitwieser et al, Analysis of Labeled Quantitative Mass Spectrometry Proteomics Data, 2012
	# p. 86
	# Carrillo et al. (2010) tested different ways to summarize data: average of ratios, 
	# Libra ratio, Linear regression on inten- sities, PCA, Ratio of Sum of Intensities,
	# and Total Least Squares. They found the error to be the smallest with the sum of intensities.
	# Carrillo B, Yanofsky C, Laboissiere S, Nadon R, Kearney RE (2010)
	#  Methods for combining peptide intensities to estimate relative protein abundance. 
	#  Bioinformatics 26(1):98Ð103
	
	if(ncol(EXPDESIGNOBJ$expDesign) > 1){ ### calc ratios if more than 1 condition
		RATIOSMEDIANNORMINTPERCOND <-  getRatiosPerCondition(MEDIANNORMINTDATA
				,controlCondition=CONTROLCONDITION
				,nonControlConditions=NONCONTROLCONDITIONS)
		
	}else{
		RATIOSMEDIANNORMINTPERCOND <- rep(NA,nrow(MEDIANNORMINTDATA))
	}
	
	
	
	### create data frame of DE qvalues
	if( ISDEEXP) { # only if each condition has at least 2 replicates or 1 condition
		
		### LOAD EXT LIBRARIES
		suppressPackageStartupMessages(library(limma,quiet=TRUE))
		suppressPackageStartupMessages(library(affy,quiet=TRUE))
		suppressPackageStartupMessages(library(gplots, quiet=TRUE ))
		suppressPackageStartupMessages(library(RSVGTipsDevice, quiet=TRUE ))
		suppressPackageStartupMessages(library(seqinr, quiet=TRUE ))
		### LOAD EXT LIBRARIES END
		
		### calculate cv per condition
		CVSPERCOND <- getCVPerCondition(GROUPEDNORMINTDATA)
	
	
	### Breitwieser et al, Analysis of Labeled Quantitative Mass Spectrometry Proteomics Data, 2012
	# p. 87
	# 5.3.6.4 Statistical Power
	# What is the number of samples required to be able to observe a certain fold change,
	# given the overall data variability? Using simulation experiments, Levin (2011) showed 
	# that when the combined technical and biological variation is as low as 25%, a fold 
	# change of 1.5 can be measured reliably with four biological replicates per sample group.
	# Levin Y (2011) The role of statistical power analysis in quantitative proteomics. Proteomics 11(12):2565Ð2567
	
	# MY COMMENT: These numbers correspond to a pval cut-off of 0.05. If we additionally apply a fold change cut-off 
	# the statistical power is of course a lot worse.
	
		statTest <- .testDE(GROUPEDNORMINTDATA
				, controlCondition=CONTROLCONDITION
				, nonControlConditions=NONCONTROLCONDITIONS
				, expDesign=EXPDESIGNOBJ$expDesign
				, pairWise=!userOptions$nonPairWiseStatTest
				, selection <- 	apply(CVSPERCOND,1,max) < userOptions$cvCutOff 
		)
		
		QVALUESPERCOND <- 		statTest$qvaluesPerCond
		PVALUESPERCOND <- 		statTest$pvaluesPerCond
		
	}else{
		cat("WARNING: NOT ENOUGH REPLICATES PER CONDITION -> NO QVALUES CALCULATED\n")
		### TURN OFF GRAPHICS
		userOptions$isVolcanoPlots <- FALSE
		#userOptions$isIntensityDistributionPlots <- FALSE
		userOptions$isHClustPlot <- FALSE
		userOptions$isDeFdrPlot <- FALSE
		userOptions$isCreateSVGFigs <- FALSE
	}
	
	###  STAT. ANALYSIS, CALCULATE RATIOS AND QVALUES OF DE FEATURES, CVS, MARK CONDITIONS CONTAINING MISSING VALUES END

###################################### EXPRESSION ANALYSIS END


###### OUTPUT #####################################################################

######################################## GRAPHICS
	### PDF OUTPUT
	source(paste(BASEDIR,"exec/subscript/outputPdfGraphics.R", sep=""))
	### PDF OUTPUT END
	### SVG OUTPUT
	source(paste(BASEDIR,"exec/subscript/outputSvgGraphics.R", sep=""))
	### SVG OUTPUT END

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
	specCountSumExport <- data.frame(getSummedIdxPerCond(getCondGroupedSampleData(UNGROUPEDSPECCOUNTDATA,expDesign=EXPDESIGNOBJ$expDesign )))
	names(specCountSumExport) <- paste("specCountSum_",names(EXPDESIGNOBJ$expDesign ),sep="")
	
	### median ratios
	medianRatiosExport <- data.frame(RATIOSMEDIANNORMINTPERCOND)
	names(medianRatiosExport) <- paste("medianRatio_",names(RATIOSMEDIANNORMINTPERCOND),"_vs_",CONTROLCONDITION,sep="")
	
	if(ISPEPTIDEANALYSIS){
		
		export <- data.frame(PROGOUTDATA$peptides
				,PROGOUTDATA$modifs
				,PROGOUTDATA$minMassDiff
				,PROGOUTDATA$Accession
				,PROGOUTDATA$Description
				,PROGOUTDATA$Confidence.score
				,PROGOUTDATA$mz
				,PROGOUTDATA$allMzs
				,PROGOUTDATA$chargeStates
				,PROGOUTDATA$retentionTime
				,rawDataExport
				,normDataExport
				,medianNormIntExport
				,specCountSumExport
				,medianRatiosExport
		)
		
		if(!is.na(userOptions$proteinFastaFile)){
			
			cat("EXTRACTING MODIFICATION SITE COORDINATES FROM ",  userOptions$proteinFastaFile,"\n")
			export <- cbind(export,getPhosphoCoordExport(PROGOUTDATA,userOptions$proteinFastaFile))
			
		}
		
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
	if(ISDEEXP){
		
		### d.e qvalues
		qvaluesExport <- data.frame(QVALUESPERCOND)
		names(qvaluesExport) <- paste("qvalue_",names(QVALUESPERCOND),"_vs_",CONTROLCONDITION,sep="")
		
		### d.e pvalues
		pvaluesExport <- data.frame(PVALUESPERCOND)
		names(pvaluesExport) <- paste("pvalue_",names(PVALUESPERCOND),"_vs_",CONTROLCONDITION,sep="")
		
		### coefficient of variance (%)
		cvsExport <- data.frame(CVSPERCOND*100)
		names(cvsExport) <- paste("cv_",names(CVSPERCOND),sep="")
		
		### experimental
		if(userOptions$eBayes){
			### add to export
			export <- data.frame(export,qvaluesExport,pvaluesExport,cvsExport)	
		}else{
			### add to export
			export <- data.frame(export,qvaluesExport,cvsExport)
		}
	}
	
	### clean up column names
	names(export) <- gsub("PROGOUTDATA.","",names(export))
	
	write.table(export,file=userOptions$csvFilePath, sep="\t", row.names=FALSE)
	
	cat("CREATED FILE ",  userOptions$csvFilePath,"\n")

###################################### TSV FILE EXPORT END


###################################### SAVE SERIALIZED R-DATA OBJECT

	if(userOptions$isSaveRObject){
		
		GROUPEDRAWINTDATA <- getCondGroupedSampleData(UNGROUPEDRAWINTDATA,expDesign=EXPDESIGNOBJ$expDesign )
		
		GROUPEDSPECCOUNTDATA <- getCondGroupedSampleData(UNGROUPEDSPECCOUNTDATA,expDesign=EXPDESIGNOBJ$expDesign )
		
		PROTEINS <- (PROGOUTDATA$Accession)
	
		if(min(EXPDESIGNOBJ$expDesign) > 1){
			save(PROTEINS
			     	,GROUPEDRAWINTDATA
					,GROUPEDSPECCOUNTDATA
					,GROUPEDNORMINTDATA
					,MEDIANNORMINTDATA
					,RATIOSMEDIANNORMINTPERCOND
					,QVALUESPERCOND
					,CVSPERCOND		
					,CONTROLCONDITION
					,EXPDESIGNOBJ 	
					,PROGOUTDATA ### TMP
					,file=userOptions$rDataFile)
		}else{ ### no qvals, cvs ... to  export
			save(PROTEINS
		     		,GROUPEDRAWINTDATA
					,GROUPEDSPECCOUNTDATA
					,GROUPEDNORMINTDATA
					,MEDIANNORMINTDATA
					,RATIOSMEDIANNORMINTPERCOND
					,CONTROLCONDITION
					,EXPDESIGNOBJ 	
					,PROGOUTDATA ### TMP
					,file=userOptions$rDataFile)
		}	
		cat("CREATED FILE ", userOptions$rDataFile,"\n")
	}

###################################### SAVE SERIALIZED R-DATA OBJECT END 

###################################### TOP X EXPORT

	if(ISPEPTIDEANALYSIS && (userOptions$topX > 0)){
		source(paste(BASEDIR,"R/AbsoluteQuantification.R", sep=""))
		
		topXDf <- calculatePQIPerProteinAndSample2(UNGROUPEDNORMINTDATA,PROGOUTDATA$Accession, method=userOptions$topX)
		topXGrouped <- getCondGroupedSampleData(topXDf,expDesign=EXPDESIGNOBJ$expDesign )
		
		write.table(file=userOptions$topXCsvFile,
				cbind(rownames(topXDf)
						,data.frame(unique(PROGOUTDATA$Description),row.names=unique(PROGOUTDATA$Accession))[rownames(topXDf),] ### add desciptions
						,topXDf
						,getMedianIdxPerCond(topXGrouped)
						,getCVPerCondition(topXGrouped)
				)
						, sep="\t",row.names=F 
				,col.names=c("Protein","Description",names(UNGROUPEDRAWINTDATA),paste("median","_",names(GROUPEDNORMINTDATA),sep=""),paste("cv","_",names(GROUPEDNORMINTDATA),sep=""))
			)
		cat("CREATED FILE ", userOptions$topXCsvFile,"\n")
	}

###################################### TOP X EXPORT END

###### OUTPUT END #####################################################################

cat("DONE\n")

