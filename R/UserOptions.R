# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### VERBOSE	
#	-v verbose 

## I/O
#	-i inputFile
#	-o outputdir
#	-l resultFilesLabel

# FILTER (--F)
#	--FP FProteinAccessionSelection
#	--FM FModificationSelection (LFQ PEP ONLY)
#   --FF FFdrCutOff  (LFQ ONLY)
#	--FD FDeltaMassTolerancePrecursor (LFQ PEP ONLY)
#	--FC FCoefficientOfVarianceMax
#  	--FN FNumberOfPeptidesPerProteinMin # protein script specific

# STATISTICS (--S)

#	--SA SAnchorProtein  
#	--SP SPvalueInclude
#		(present Moderated t-statistic p-values from eBayes)

# EXPERIMENTAL DESIGN (--E)

#	--EX EXperimentalDesign

# PDF-REPORT (--P)
#	--PR PRatioCutOff
#	--PQ PQvalueCutOff
#	--PS PSelectedGraphics

# TSV-REPORT (--T)
#	--TF	TFastaFile
#	--TP	TProtease

# ADDITIONAL-REPORTS (--A)
#	--AR	ARDataFile
#	--AI	AIbaq
#	--AT	ATopX


## TEST
# -t


### LOAD EXT LIBRARIES
suppressPackageStartupMessages(library(optparse))
### LOAD EXT LIBRARIES END
### CMD OPTIONS
option_list <- list(
		
		
### I/O
		make_option(c("-i", "--inputFile"), type="character", default="",
				help="I/O: Progenesis csv input file (REQUIRED)",
		),	
		make_option(c("-o", "--outputDir"), type="character", default="",
				help="I/O:  Results Output Directory [default ./]",
		),
		
		make_option(c("-l", "--resultsFileLabel"), type="character", default="analysisResults",
				help="I/O: .pdf & .tsv file labels (prefix) [default %default]", 
		),

### I/O END
		
# FILTER (--F)
		make_option(c("--FProteinAccessionSelection"), type="character", default=".",
				help="FILTER: --FP Filter features by Accession Regular Expression [default %default] (all features kept)",
				metavar="Protein Accession Reg. expr."),
		
		#### peptide script specfic
		make_option(c("--FModificationSelection"), type="character", default=".",
				help="FILTER (LFQ PEP ONLY): --FM Only keep Peptides with modifications matching Regular Expression [default %default]
				 (all features kept). Peptide script ONLY",
				metavar="modification name Reg. expr."),
		
		
		make_option(c("--FFdrCutoff"), type="double", default=0.01,
				help="FILTER (LFQ ONLY): --FF Identification level False Discovery Rate Cutoff.  [0-1] [default %default]",
				metavar="Peptide/Protein FDR cutoff"),
		
		make_option(c("--FCoefficientOfVarianceMax"), type="double", default=Inf,
				help="FILTER: --FC Do not include features with C.V. above this threshold in statistical 
				test for differential expression [default %default]",
				metavar="Coefficent of Variance cutoff"),
		
		#### peptide script specfic
		make_option(c("--FDeltaMassTolerancePrecursor"), type="character", default="[-10,10]",
				help="FILTER (LFQ PEP ONLY): --FD Precursor mass Error Range filter [default %default ppm].
				Peptide script ONLY",
				metavar="Mass Range [x,y]"),
		
		#### protein script specfic
		make_option(c("--FNumberOfPeptidesPerProteinMin"), type="integer", default=1,
				help="FILTER: --FN Only include those proteins with at least x identified peptides [default %default]
				Protein script ONLY.",
				metavar="Number of peptides"),
		
# FILTER (--F) END	
		
# STATISTICS (--S)
	
	make_option(c("--SAnchorProtein"), type="character", default=".",
			help="STATISTICS: --SA Normalize Intensities by selected protein(s) Regular Expression
			 [default %default] (use all proteins). 
			!!! Note that the specified proteins will be excluded from results output !!!",
			metavar="Protein Accession Reg. expr."),
	
	make_option(c("--SPvalueInclude"), action="store_true", default=FALSE,
			help="STATISTICS: --SP output eBayes moderated t-statistic p-values [default %default]"),
# STATISTICS (--S) END

# EXPERIMENTAL DESIGN (--E)

	make_option(c("--EXperimentalDesign"), type="character", default=NA,
			help='EXPERIMENTAL DESIGN: --EX "," seperated samples, ":" separated conditions 
					Example: 1,2,3:4,5,6 
					->  condition1 (REF) : channel 1,2,3
					->  condition2: channel 4,5,6
					Note: for 10-plex default is "1,4,7,10:2,5,8:3,6,9"
					[default %default]'), 
# EXPERIMENTAL DESIGN (--E) END

# PDF-REPORT (--P) 
	make_option(c("--PRatioCutOff"), type="double", default=1,
		help="PDF-REPORT: --PR Intensity ratio (fold change) cut-off used for graphics export. >1 [default %default]",
		metavar="Intensity ratio cutoff"),	

	make_option(c("--PQvaueCutOff"), type="double", default=0.01,
			help="PDF-REPORT: --PQ Qvalue cut-off used for graphics. 
			High-lighting features with a qval < specified value. [0-1] [default %default]",
			metavar="Differential expression qvalue cutOff"),	
	
	make_option(c("--PSelectedGraphics"), type="character", default="",
			help="PDF-REPORT: --PS Excluded Graphics: give letter for each plot to exclude ex: --PS iv 
					(creates all plots but intensity density plots & volcano plot)
					experimental design (e)
					peptide feature score distrib related plots (f)
					intensity distibution plots (i)
					volcano plots (v)
					hierarchical clustering plots (h)
					differential expression fdr plot (d)	
					[default (all plots) %default]"),		
# PDF-REPORT (--P) END

# TSV-REPORT (--T)
	
	make_option(c("--TFastaFile"), type="character", default="",
			help="TSV-REPORT (LFQ PEP): -TF Protein Fasta File used to extract Modification Site Coordinates [default None]",
			metavar=".fasta file path"),	
# TSV-REPORT (--T) END

# ADDITIONAL-REPORTS (--A)
	make_option(c("--ARDataFile"), action="store_true", default=FALSE,
		help="ADDITIONAL-REPORTS: --AR Save R objects in 'label'.RData file [default %default]"),

	make_option(c("--AProtease"), type="character", default="KR",
			help="ADDITIONAL-REPORTS: --TP protease [default (trypsin) %default]
					1) trypsin
					2) lys-c
			
					Option considered for iBAQ normalization", 
	),

	make_option(c("--AIbaq"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS (LFQ PROT): --AI creates .tsv output file
					including protein iBAQ values. [default %default]"),
	
	make_option(c("--ATop3"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS (LFQ PEP): --AT creates .tsv output file
					including protein top3 values. [default %default]"),
	
# ADDITIONAL-REPORTS (--A) END

# TEST (peptide script specific)
	make_option(c("-t", "--test"), action="store_true", default=FALSE,
			help="TEST: test option, include first 2000 entries only [default %default]
			Peptide script ONLY."),
# TEST END
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
			help="Print extra output [default %default]")
	)

	
#' Read User Specified Command Line Options
#' @param version Safequant version number
#' @return user options list
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getUserOptions <- function(version=version){
	
	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser( prog=paste("SafeQuant",version), option_list=option_list))
	
	### CMD OPTIONS END						
	
	### SET USER OPTIONS
	userOptions <- list()

### VERBOSE
	#VERBOSE: verbose
	userOptions$verbose <- cmdOpt$verbose
### VERBOSE	END
	
# I/O
	#I/O: progenesisFilePath
	userOptions$inputFile <- cmdOpt$inputFile
	if( userOptions$inputFile == "" | !file.exists(userOptions$inputFile)){
		cat("ERROR. Please specify input file.",userOptions$inputFile, "Not found!","\n")
		q(status=-1)
	}
	
	#I/O: outputDir
	userOptions$outputDir <- cmdOpt$outputDir
	if(!file.exists(userOptions$outputDir) & userOptions$outputDir != "" ){
		cat("ERROR. No such directory",userOptions$outputDir,"\n")
		q(status=-1)
	}else if(substr(userOptions$outputDir,nchar(userOptions$outputDir),nchar(userOptions$outputDir)) != "/"){
		if(userOptions$verbose){
			cat("added slash to outputdir","\n")
		}
		if(!(userOptions$outputDir== "")){
			userOptions$outputDir <- paste(userOptions$outputDir,"/", sep="")
		}
	}
	
	#I/O: resultsFileLabel
	userOptions$resultsFileLabel <- cmdOpt$resultsFileLabel
	
	#I/O: set export file paths
	userOptions$pdfFilePath <- paste(userOptions$outputDir,userOptions$resultsFileLabel,".pdf",sep="")
	userOptions$tsvFilePath <- paste(userOptions$outputDir,userOptions$resultsFileLabel,".tsv",sep="")
	
# I/O END
	
# FILTER (--F)

	#FILTER: selectedProteinName
	userOptions$selectedProteinName <- cmdOpt$FProteinAccessionSelection
	
	#FILTER: selectedModifName
	userOptions$selectedModifName <- cmdOpt$FModificationSelection
	
	#FILTER: fdrCutoff
	userOptions$fdrCutoff <- cmdOpt$FFdrCutoff
	if(is.na(userOptions$fdrCutoff) | userOptions$fdrCutoff <= 0 | userOptions$fdrCutoff > 1 ){
		cat("ERROR. fdrCutOff must be in the range [0-1]. You specified",userOptions$fdrCutoff,"\n")
		q(status=-1)
	}
	
	#FILTER: precursorMassFilter
	### check input format precursorMassFilter
	if(is.na(cmdOpt$FDeltaMassTolerancePrecursor) | regexpr("^\\[\\-*[0-9]{1,}\\,\\-*[0-9]{1,}\\]$",as.character(cmdOpt$FDeltaMassTolerancePrecursor)) == -1 ){
		cat("ERROR. invalid precursorMassFilter ",userOptions$precursorMassFilter,"\n") 
		q(status=-1)
	}
	### add lower and upper mass bound to vector
	userOptions$precursorMassFilter <- gsub("(\\[)","",cmdOpt$FDeltaMassTolerancePrecursor)
	userOptions$precursorMassFilter <- gsub("(\\])","",userOptions$precursorMassFilter)
	userOptions$precursorMassFilter <- sort(as.numeric(unlist(strsplit(userOptions$precursorMassFilter,","))))
	
	#FILTER: cvCutOff
	userOptions$cvCutOff <- cmdOpt$FCoefficientOfVarianceMax
	if(is.na(userOptions$cvCutOff) | (userOptions$cvCutOff < 0)){
		print(paste("ERROR. cvCutOff must be > 0. You specified ", userOptions$cvCutOff))
		q(status=-1)
	}

	#FILTER: minNbPeptidesPerProt
	userOptions$minNbPeptidesPerProt <- cmdOpt$FNumberOfPeptidesPerProteinMin
	if(is.na(userOptions$minNbPeptidesPerProt) | userOptions$minNbPeptidesPerProt < 0 ){
		print(paste("ERROR. peptidesForQuantCutoff must be >= 0. You specified ", userOptions$minNbPeptidesPerProt))
		q(status=-1)
	}
	
# FILTER (--F) END

# STATISTICS
	
	#STATISTICS: normAC
	userOptions$normAC <- cmdOpt$SAnchorProtein

	#STATISTICS: eBayes
	userOptions$eBayes <- cmdOpt$SPvalueInclude
	
# STATISTICS END	
	
# EXPERIMENTAL DESIGN

	#EXPERIMENTAL DESIGN: EXperimentalDesign
	userOptions$expDesignTag <- cmdOpt$EXperimentalDesign
	
# EXPERIMENTAL DESIGN END

# PDF-REPORT (--P)

	# PDF-REPORT: ratioCutOff
	userOptions$ratioCutOff <- cmdOpt$PRatioCutOff
	if(is.na(userOptions$ratioCutOff) | userOptions$ratioCutOff < 1){
		cat("ERROR. ratioCutoff must be > 1. You specified",userOptions$ratioCutOff,"\n")
		q(status=-1)
	}
	
	# PDF-REPORT: deFdrCutoff
	userOptions$deFdrCutoff <- cmdOpt$PQvaueCutOff
	if(is.na(userOptions$deFdrCutoff) | userOptions$deFdrCutoff <= 0 | userOptions$deFdrCutoff > 1 ){
		cat("ERROR. deFdrCutoff must be in the range [0-1]. You specified",userOptions$deFdrCutoff,"\n")
		q(status=-1)
	}

	# PDF-REPORT: PSelectedGraphics
	userOptions$isDispExpDesign <- !regexpr("e",cmdOpt$PSelectedGraphics) > -1
	userOptions$isFdrPlots <- !regexpr("f",cmdOpt$PSelectedGraphics) > -1
	userOptions$isIntensityDistributionPlots <- !regexpr("i",cmdOpt$PSelectedGraphics) > -1
	userOptions$isVolcanoPlots <- !regexpr("v",cmdOpt$PSelectedGraphics) > -1
	userOptions$isHClustPlot <- !regexpr("h",cmdOpt$PSelectedGraphics) > -1
	userOptions$isDeFdrPlot <- !regexpr("d",cmdOpt$PSelectedGraphics) > -1


# PDF-REPORT (--P) END

# TSV-REPORT (--T)
	
	# TSV-REPORT: proteinFastaFile
	userOptions$proteinFastaFile <- NA
	if(nchar(cmdOpt$TFastaFile) > 0 ){
		### check if file exists
		if(file.exists(cmdOpt$TFastaFile)){
			userOptions$proteinFastaFile <- cmdOpt$TFastaFile
		}else{
			cat("ERROR. File does not exist",cmdOpt$TFastaFile,"\n")
			q(status=-1)
		}				
	}

	# TSV-REPORT: proteaseTarget	(deprecated)
	userOptions$protease <- cmdOpt$TProtease	

# TSV-REPORT (--T) END

# ADDITIONAL-REPORTS (--A)

	#ADDITIONAL-REPORTS iBaqTsvFile
	userOptions$iBAQ <- cmdOpt$AIbaq
	userOptions$iBAQFile <- paste(userOptions$outputDir,userOptions$resultsFileLabel,"_iBAQ.tsv",sep="")

    #ADDITIONAL-REPORTS top3TsvFile
	userOptions$top3 <- cmdOpt$ATop3
	userOptions$top3File <- paste(userOptions$outputDir,userOptions$resultsFileLabel,"_top3.tsv",sep="")

	#ADDITIONAL-REPORTS rDataFile, isSaveRObject	
	userOptions$isSaveRObject <- cmdOpt$ARDataFile
	if(userOptions$isSaveRObject){
		userOptions$rDataFile <- paste(userOptions$outputDir,userOptions$resultsFileLabel,".rData",sep="")
	}
	
	
# ADDITIONAL-REPORTS (--A) END


# TEST

	### test run to define parameters (peptide script specific)
	userOptions$test <- cmdOpt$test

# TEST END
	return(userOptions)

}


userInputTag <- "1,2,3:4,5,6"

# tag: 1,2:3:4,5,6 
#condition isControl
#1 Condition 1      TRUE
#2 Condition 1      TRUE
#3 Condition 1     TRUE
#4 Condition 2     FALSE
#5 Condition 2     FALSE
#6 Condition 2     FALSE

#' Create experimental design data.frame from user input string
#' @param string tag
#' @param data.frame default expDesign 
#' @return data.frame describing experimental design
#' @export
#' @note  No note
#' @details  tag: 1,2:3:4,5,6 
#'		condition isControl
#'	1 Condition 1 TRUE
#'	2 Condition 1 TRUE
#'	3 Condition 1 TRUE
#'	4 Condition 2 FALSE
#'	5 Condition 2 FALSE
#'	6 Condition 2 FALSE
#' @references NA 
#' @examples print("No examples")
expDesignTagToExpDesign <- function(tag, expDesignDefault){
	sampleOrder <- as.numeric(unlist(strsplit(tag,"[\\,\\:]")))

	# make sure no duplicates, withing range etc.
	if(is.na(sampleOrder[1]) 
			| (max(table(sampleOrder))>1) 
			| (min(sampleOrder) < 1)
			| max(as.numeric(sampleOrder)) > nrow(expDesignDefault)
			){
		stop("ERROR: getExpDesign, INVALID EXPERIMENTAL DESIGN ",tag,"\n")
		
	}
	
	expDesign <- data.frame(row.names=sampleOrder,condition=rep(NA,length(sampleOrder)), isControl=rep(FALSE,length(sampleOrder))  )
	
	condNb <- 1
	for(cond in unlist(strsplit(tag,":"))){
		
		expDesign[as.character(unlist(strsplit(cond,","))),]$condition <- paste("Condition",condNb)
		#expDesign <- cbind(expDesign,length(unlist(strsplit(cond,","))))
		
		condNb <- condNb + 1
	}       
	
	expDesign[ expDesign[,1] == "Condition 1" ,]$isControl <- T
	
	### get original sample names
	rownames(expDesign) <- rownames(expDesignDefault)[as.numeric(rownames(expDesign))]
	#expDesignUser$condition <- expDesign[rownames(expDesignUser) ,]$condition
	
	return(expDesign)
	
}


