#!/usr/bin/Rscript

# TODO: Add comment
# 
# Author: erikahrne
###############################################################################
#### DEPENDANCIES
#limma
#affy
#gplots
#optparse
#RSVGTipsDevice


### GET TMTDIR & BASEDIR
initial.options <- commandArgs(trailingOnly = FALSE)
BASEDIR <- paste(normalizePath(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))),"/",sep="")
BASEDIR <- gsub("TMT","",BASEDIR,)
BASEDIR <- gsub("exec","",BASEDIR)
BASEDIR <- gsub("\\/{2,2}","",BASEDIR)
BASEDIR <- gsub("\\\\{2,2}","",BASEDIR)

remove(initial.options)
### GET TMTDIR & BASEDIR END

### SOURCE GENERAL FUNCTIONS
#### LOAD IDENTIFICATION RELATED FUNCTIONS (used by CsvParser)
source(paste(BASEDIR,"R/IdentificationAnalysis.R", sep=""))
### LOAD EXPRESSION ANALYSIS FUNCTIONS
source(paste(BASEDIR,"R/ExpressionAnalysis.R", sep=""))
### LOAD PLOTTING FUNCTIONS
source(paste(BASEDIR,"R/Graphics.R", sep=""))
### LOAD SVG PLOTTING FUNCTIONS
source(paste(BASEDIR,"R/SVGPlots.R", sep=""))

### LOAD TMT FUNCTIONS
source(paste(BASEDIR,"R/TMT.R", sep=""))

### GET VERSION
source(paste(BASEDIR,"config/VERSION.R", sep=""))
### END GET VERSION


### USER CMD LINE OPTIONS
source(paste(BASEDIR,"R/userOptions.R", sep=""))
userOptions <- getUserOptions(version=VERSION)

### SUPRESS WARNINGS
if(!userOptions$verbose){
	options(warn=-1)
}

### PARSE SCAFFOLD RAW DATA FILE
	SCAFFOLDOUTDATA <- parseScaffoldRawFile(userOptions$inputFilePath)
	### ac -> protein descrption dictionary
	PROTEINDESCDIC <- unique(data.frame(SCAFFOLDOUTDATA$Protein.Name,SCAFFOLDOUTDATA$Accession.Numbers))
	row.names(PROTEINDESCDIC) <- as.character(PROTEINDESCDIC[,2])
	
	### test mode, sample 300 feature entries 
	if(userOptions$test){
		cat("WARNING: SafeQuant is run in TEST mode \n")
		SCAFFOLDOUTDATA <- SCAFFOLDOUTDATA[sample(1:nrow(SCAFFOLDOUTDATA),300),]
	}
	
### PARSE SCAFFOLD RAW DATA FILE END

### GET EXPERIMENT DESIGN

	NBPLEX <- getNbPlex(SCAFFOLDOUTDATA)

	# if 10-plex (and default exp tag) 	
	if(NBPLEX == 10 && userOptions$expDesignTag == '1,3,5:2,4,6' ){
		userOptions$expDesignTag <- "1,4,7,10:2,5,8:3,6,9"
		if(userOptions$verbose){
			cat("INFO: 10-plex Experiment \n")
		}
	}
	
	EXPDESIGNOBJ <- getExpDesignObj(userOptions$expDesignTag,NBPLEX)

	### CHECK IF MORE THAN ONE CONDITION WITH MORE THAN ONE REP (IF NOT NO STAT ANALYSIS AND ACCOMPANYING POTS)
	ISDEEXP <- (min(EXPDESIGNOBJ$expDesign) > 1) & (ncol(EXPDESIGNOBJ$expDesign) > 1) ### @TODO get rid off this, add test to each function
	
	if(userOptions$verbose){
		cat("INFO expDesign\n")
		print(EXPDESIGNOBJ)
	}
	
### GET EXPERIMENT DESIGN END

### FILTER

	# filter by protein AC FProteinAccessionSelection
	SCAFFOLDOUTDATA <- SCAFFOLDOUTDATA[regexpr(userOptions$selectedProteinName,SCAFFOLDOUTDATA$Accession.Numbers,ignore.case = TRUE ) > -1 ,]

### FILTER END

######################################## EXPRESSION ANALYSIS

	

	UNGROUPEDRAWINTDATA <- SCAFFOLDOUTDATA[,10:(9+NBPLEX)]
	### channel impurity correct
	#UNGROUPEDRAWINTDATA <- data.frame(2^purityCorrectTMT(log2(as.matrix(UNGROUPEDRAWINTDATA)),impurityMatrix=getImpuritiesMatrix()))
	UNGROUPEDRAWINTDATA <- data.frame(purityCorrectTMT(as.matrix(UNGROUPEDRAWINTDATA),impurityMatrix=getImpuritiesMatrix(NBPLEX)))
	names(UNGROUPEDRAWINTDATA) <- names(SCAFFOLDOUTDATA[,10:(9+NBPLEX)])
	
	### replace missing values
	BASELINEINTENSITY <- getBaselineIntensity(unlist(UNGROUPEDRAWINTDATA))
	UNGROUPEDRAWINTDATA[UNGROUPEDRAWINTDATA == 0] <- BASELINEINTENSITY
	
	# protein level
	proteinSummary <- getIntSumPerProtein(UNGROUPEDRAWINTDATA,SCAFFOLDOUTDATA$Accession.Numbers,SCAFFOLDOUTDATA$Peptide.Sequence,minNbPeptPerProt= userOptions$minNbPeptidesPerProt)
	UNGROUPEDRAWINTDATA <- proteinSummary$perProteinIntSum
	
	# update sample order after purity correction	
	UNGROUPEDRAWINTDATA <- data.frame(UNGROUPEDRAWINTDATA[,EXPDESIGNOBJ$sampleOrder])
	names(UNGROUPEDRAWINTDATA) <- names(proteinSummary$perProteinIntSum)[EXPDESIGNOBJ$sampleOrder] ## in case just one sample selected
	
	# normalize	
	normFactors <- getNormalizationFactors(UNGROUPEDRAWINTDATA[( regexpr(userOptions$normAC,rownames(UNGROUPEDRAWINTDATA),ignore.case=TRUE) > -1),]) ### allow regexpr of user selecte norm proteins
	UNGROUPEDNORMINTDATA <- normalizeIntensities(UNGROUPEDRAWINTDATA, normFactors)
	
	GROUPEDNORMINTDATA <- getCondGroupedSampleData(UNGROUPEDNORMINTDATA,expDesign=EXPDESIGNOBJ$expDesign )
	names(GROUPEDNORMINTDATA) <- names(EXPDESIGNOBJ$expDesign)
	
	CONTROLCONDITION <- names(EXPDESIGNOBJ$expDesign)[1]
	NONCONTROLCONDITIONS <- names(EXPDESIGNOBJ$expDesign)[2:ncol(EXPDESIGNOBJ$expDesign)]
	
	### create data frame of ratios of medians
	MEDIANNORMINTDATA <- getMedianIdxPerCond(GROUPEDNORMINTDATA)
	
	### calculate ratios
	if(ncol(EXPDESIGNOBJ$expDesign) > 1){ ### calc ratios if more than 1 condition
		RATIOSMEDIANNORMINTPERCOND <-  getRatiosPerCondition(MEDIANNORMINTDATA
				,controlCondition=CONTROLCONDITION
				,nonControlConditions=NONCONTROLCONDITIONS)
		
	}else{
		RATIOSMEDIANNORMINTPERCOND <- rep(NA,nrow(MEDIANNORMINTDATA))
	}
	
	if(ISDEEXP){
		
		CVSPERCOND <- getCVPerCondition(GROUPEDNORMINTDATA)
		
		### LOAD EXT LIBRARIES
		suppressPackageStartupMessages(library(limma,quiet=TRUE))
		suppressPackageStartupMessages(library(affy,quiet=TRUE))
		suppressPackageStartupMessages(library(gplots, quiet=TRUE ))
		suppressPackageStartupMessages(library(RSVGTipsDevice, quiet=TRUE ))
		suppressPackageStartupMessages(library(seqinr, quiet=TRUE ))
		### LOAD EXT LIBRARIES END
		
		statTestSummed <- .testDE(GROUPEDNORMINTDATA
				, controlCondition=CONTROLCONDITION
				, nonControlConditions=NONCONTROLCONDITIONS
				, expDesign=EXPDESIGNOBJ$expDesign
				, pairWise=T
				, selection <- 	apply(CVSPERCOND,1,max) < userOptions$cvCutOff 
		)
		
		QVALUESPERCOND <- 		statTestSummed$qvaluesPerCond
		PVALUESPERCOND <- 		statTestSummed$pvaluesPerCond
		
	}else{
		cat("WARNING: NOT ENOUGH REPLICATES PER CONDITION -> NO QVALUES CALCULATED\n")
		### TURN OFF GRAPHICS
		userOptions$isVolcanoPlots <- FALSE
		#userOptions$isIntensityDistributionPlots <- FALSE
		userOptions$isHClustPlot <- FALSE
		userOptions$isDeFdrPlot <- FALSE
		userOptions$isCreateSVGFigs <- FALSE
	
	}
	


######################################## EXPRESSION ANALYSIS END


###### OUTPUT #####################################################################

###################################### TSV FILE EXPORT

### TSV EXPORT
### @TODO 
	
export <- data.frame(rownames(UNGROUPEDNORMINTDATA)
		,PROTEINDESCDIC[rownames(UNGROUPEDNORMINTDATA),1]
		,proteinSummary$peptidesPerProtein	
		,proteinSummary$spectraPerProtein	
		,MEDIANNORMINTDATA
		,UNGROUPEDNORMINTDATA
		,log2(RATIOSMEDIANNORMINTPERCOND)

)

names(export)[1:4] <- c("AC","description","Peptides_Per_Protein","SSM_Per_Protein")				
names(export)[5:(5+EXPDESIGNOBJ$nbConditions-1)] <- paste("medianNormInt_cond_",1:EXPDESIGNOBJ$nbConditions,sep="") 
i <- 5+EXPDESIGNOBJ$nbConditions
#names(export)[i:(i+EXPDESIGNOBJ$nbConditions-1)] <- paste("cv_cond",1:EXPDESIGNOBJ$nbConditions,sep="")
#i <- i+EXPDESIGNOBJ$nbConditions
names(export)[i:ncol(export)] <- gsub("X","Channel_",names(export)[i:ncol(export)]) 
i <- i+sum(EXPDESIGNOBJ$expDesign)
names(export)[i:(i+EXPDESIGNOBJ$nbConditions-2)] <- paste("log2_ratio_cond_",((1:(EXPDESIGNOBJ$nbConditions-1))+1),"_vs_",CONTROLCONDITION,sep="") 

### qvalues and coef. variance not calculated if only one sample per condition 
if(ISDEEXP){
	
	### d.e qvalues
	qvaluesExport <- data.frame(QVALUESPERCOND)
	names(qvaluesExport) <- paste("qvalue_cond_",names(QVALUESPERCOND),"_vs_cond_",CONTROLCONDITION,sep="")
	
	### d.e pvalues
	pvaluesExport <- data.frame(PVALUESPERCOND)
	names(pvaluesExport) <- paste("pvalue_cond_",names(PVALUESPERCOND),"_vs_cond_",CONTROLCONDITION,sep="")
	
	### coefficient of variance (%)
	cvsExport <- data.frame(CVSPERCOND*100)
	names(cvsExport) <- paste("cv_cond_",names(CVSPERCOND),sep="")
	
	### experimental
	if(userOptions$eBayes){
		### add to export
		export <- data.frame(export,qvaluesExport,pvaluesExport,cvsExport)	
	}else{
		### add to export
		export <- data.frame(export,qvaluesExport,cvsExport)
	}
}

write.table(file=userOptions$csvFilePath,export,row.names=F,sep="\t")

cat("CREATED FILE ",  userOptions$csvFilePath,"\n")

###################################### TSV FILE EXPORT END

######################################## GRAPHICS
### PDF OUTPUT
source(paste(BASEDIR,"exec/subscript/outputPdfGraphics.R", sep=""))
### PDF OUTPUT END
### SVG OUTPUT
source(paste(BASEDIR,"exec/subscript/outputSvgGraphics.R", sep=""))
### SVG OUTPUT END

######################################## GRAPHICS END

###################################### SAVE SERIALIZED R-DATA OBJECT

if(userOptions$isSaveRObject){
	
	GROUPEDRAWINTDATA <- getCondGroupedSampleData(UNGROUPEDRAWINTDATA,expDesign=EXPDESIGNOBJ$expDesign )
	
	if(min(EXPDESIGNOBJ$expDesign) > 1){
		save(  	GROUPEDRAWINTDATA
				,GROUPEDNORMINTDATA
				,MEDIANNORMINTDATA
				,RATIOSMEDIANNORMINTPERCOND
				,QVALUESPERCOND
				,CVSPERCOND		
				,CONTROLCONDITION
				,EXPDESIGNOBJ 	
				,SCAFFOLDOUTDATA ### TMP
				,proteinSummary
				,file=userOptions$rDataFile)
	}else{ ### no qvals, cvs ... to  export
		save(	GROUPEDRAWINTDATA
				,GROUPEDNORMINTDATA
				,MEDIANNORMINTDATA
				,RATIOSMEDIANNORMINTPERCOND
				,CONTROLCONDITION
				,EXPDESIGNOBJ 	
				,SCAFFOLDOUTDATA ### TMP
				,proteinSummary
				,file=userOptions$rDataFile)
	}	
	cat("CREATED FILE ", userOptions$rDataFile,"\n")
}

###################################### SAVE SERIALIZED R-DATA OBJECT END

###### OUTPUT END #####################################################################

