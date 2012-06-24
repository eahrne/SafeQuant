### LOAD EXT LIBRARIES
suppressPackageStartupMessages(library(optparse))
### LOAD EXT LIBRARIES END
### CMD OPTIONS

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
						a identification score corresponding to a fdr cut-off < specified value. [0-1]
						[default %default]",
				metavar="Protein level FDR cutoff"),
		
		make_option(c("-r", "--ratioCutoff"), type="double", default=1.5,
				help="Intensity ratio (fold change) cut-off used for
						graphics export. >1 [default %default]",
				metavar="Intensity ratio cutoff"),
		
		make_option(c("-q", "--deFdrCutoff"), type="double", default=0.01,
				help="fdr cut-off used for graphics. 
						High-lighting features with a qval < specified value. [0-1] [default %default]",
				metavar="Differential expression fdr cut-off"),
#		make_option(c("-b", "--logOddsCutoff"), type="double", default=2.2,
#				help="log odds cut-off used for graphics.
#						High-lighting features with a log odds ratio > specified value. [default %default]",
#				metavar="Intensity ratio logOdds cutoff"),
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
						hierarchical clustering plots (h)
						differential expression fdr plot (d)	
						[default (all plots) %default]"),
		
		#### peptide script specfic
		make_option(c("-t", "--modificationTypeFilter"), type="character", default=".",
				help="Filter features by Modification Type [default %default] (all features kept)",
				metavar="modification name regexpr"),
		make_option(c("-m", "--precursorMassFilter"), type="double", default=10,
				help="Precursor mass filter 
						min 1 ppm [default %default ppm]"),
		make_option(c("-b", "--topX"), type="integer", default=0,
				help="creates csv output including the top X 
						most intense features per protein [default %default = off]",
				metavar="top X peptides per protein"),
		make_option(c("-d", "--proteaseTarget"), type="character", default="KR",
				help=" protease target aa's [default (Trypsin) %default]", 
		),
		make_option(c("-y", "--test"), action="store_true", default=FALSE,
				help="TEST (parse only first 1000 entries)  [default %default]")

)

getUserOptions <- function(version=version){

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser( prog=paste("SafeQuant",VERSION), option_list=option_list))
	
	### CMD OPTIONS END
	
	### SET USER OPTIONS
	userOptions <- list()
	
	#@TODO update variable names in peptidesSQAnalysis.R
	
	#minIntensity
	#userOptions$minIntensity
	
	#mvaPlots              
	
	#normAC
	#userOptions$normAC
	
	#outputDir
	#userOptions$outputDirPath
	
	#pdfFile
	#userOptions$pdfFilePath
	
	#peptidesForQuantCutoff
	#userOptions$minNbPeptidesPerProt
	
	#progenesisFile
	#userOptions$progenesisFilePath
	
	#ratioCutOff           
	#userOptions$ratioCutOff
	
	#resultsFileLabel
	#userOptions$resultsFileLabel
	
	#saveRObject
	#userOptions$isSaveRObject
	
	#selectedControlCond
	#userOptions$selectedControlCond
	
	#selectedGraphics
	
	#selectedProteinName
	#userOptions$selectedProteinName
	
	#verbose
	#userOptions$verbose
	
	#volcanoPlots
	#userOptions$isVolcanoPlots
	
	#topX
	#userOptions$topX
	#
	#selectedModifName
	#userOptions$selectedModifName
	#
	#precursorMassFilter
	#userOptions$precursorMassFilter
	
	#proteaseTarget
	#userOptions$proteaseTarget
	
	#test
	#userOptions$test
	
	### TEST REQUIRED OPTIONS
	
	userOptions$progenesisFilePath <- cmdOpt$inputFile
	if( userOptions$progenesisFilePath == "" | !file.exists(userOptions$progenesisFilePath)){
		print(paste("ERROR. Please specify input file.",userOptions$progenesisFilePath, "Not found!"))
		q(status=-1)
	}
	
	### TEST REQUIRED OPTIONS END
	
	### TEST NON-REQUIRED OPTIONS
	
	userOptions$verbose <- cmdOpt$verbose
	
	# output
	userOptions$outputDirPath <- cmdOpt$outputDir
	if(!file.exists(userOptions$outputDirPath) & userOptions$outputDirPath != "" ){
		print(paste("ERROR. No such directory",userOptions$outputDirPath))
		q(status=-1)
	}else if(substr(userOptions$outputDirPath,nchar(userOptions$outputDirPath),nchar(userOptions$outputDirPath)) != "/"){
		if(userOptions$verbose){
			print(paste("added slash to outputdir")) 
		}
		if(!(userOptions$outputDirPath== "")){
			userOptions$outputDirPath <- paste(userOptions$outputDirPath,"/", sep="")
		}
	}
	
	userOptions$resultsFileLabel <- cmdOpt$resultsFileLabel
	
	# set export file paths
	userOptions$pdfFilePath <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,".pdf",sep="")
	userOptions$csvFilePath <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,".csv",sep="")
	
	userOptions$isSaveRObject <- cmdOpt$saveRObject
	if(userOptions$isSaveRObject){
		rDataFile <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,".rData",sep="")
	}
	
	userOptions$topX <- cmdOpt$topX
	if(is.na(cmdOpt$topX) | userOptions$topX < 0){
		print(paste("ERROR. topX must be >= 0. You specified",userOptions$topX)) 
		q(status=-1)
	}
	userOptions$topXCsvFile <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,"_top",userOptions$topX,".csv",sep="")
	
	# output end
	
	# filtering
	userOptions$selectedProteinName <- cmdOpt$proteinAccesionFilter
	userOptions$normAC <- cmdOpt$proteinAccessionNormalization
	
	userOptions$minNbPeptidesPerProt <- cmdOpt$peptidesForQuantCutoff
	if(is.na(cmdOpt$peptidesForQuantCutoff) | userOptions$minNbPeptidesPerProt < 0 ){
		print(paste("ERROR. peptidesForQuantCutoff must be >= 0. You specified ", userOptions$minNbPeptidesPerProt))
		q(status=-1)
	}
	
	userOptions$fdrCutoff <- cmdOpt$fdrCutoff
	if(is.na(cmdOpt$fdrCutoff) | userOptions$fdrCutoff <= 0 | userOptions$fdrCutoff > 1 ){
		print(paste("ERROR. fdrCutOff must be in the range [0-1]. You specified",userOptions$fdrCutoff)) 
		q(status=-1)
	}
	
	userOptions$selectedModifName <- cmdOpt$modificationTypeFilter
	
	userOptions$precursorMassFilter <- cmdOpt$precursorMassFilter
	if(is.na(cmdOpt$precursorMassFilter) | userOptions$precursorMassFilter <= 0 ){
		print(paste("ERROR. precursorMassFilter must be > 1 ppm. You specified",userOptions$precursorMassFilter)) 
		q(status=-1)
	}
	
	# filtering end
	
	# data analysis
	
	### @TODO CONSIDER DISCARDING THIS OPTION (ALT. PROVIDE LIST OF COMMON PROTEASES).
	userOptions$proteaseTarget <- cmdOpt$proteaseTarget
	
	userOptions$minIntensity <- cmdOpt$naReplace
	if(is.na(cmdOpt$naReplace) | userOptions$minIntensity <= 0){
		print(paste("ERROR. naReplace must be > 0. You specified",userOptions$minIntensity)) 
		q(status=-1)
	}
	
	userOptions$selectedControlCond <- cmdOpt$controlCondition
	if(is.na(cmdOpt$controlCondition) | userOptions$selectedControlCond < 1 ){
		print(paste("ERROR. controlCondition must be >= 1. You specified ", userOptions$selectedControlCond))
		q(status=-1)
	}
	
	# data analysis end 
	
	# graphics
	
	userOptions$ratioCutOff <- cmdOpt$ratioCutoff
	if(is.na(cmdOpt$ratioCutoff) | userOptions$ratioCutOff <= 1){
		print(paste("ERROR. ratioCutoff must be > 1. You specified",userOptions$ratioCutOff)) 
		q(status=-1)
	}
	
	#BCutoff <- cmdOpt$logOddsCutoff
	#if(is.na(cmdOpt$logOddsCutoff)){
	#	print(paste("ERROR. Invalid log Odds Cutoff. You specified",BCutoff)) 
	#	q(status=-1)
	#}
	
	userOptions$deFdrCutoff <- cmdOpt$deFdrCutoff
	if(is.na(cmdOpt$deFdrCutoff) | userOptions$deFdrCutoff <= 0 | userOptions$deFdrCutoff > 1 ){
		print(paste("ERROR. deFdrCutoff must be in the range [0-1]. You specified",userOptions$deFdrCutoff)) 
		q(status=-1)
	}
	
	#selectedGraphics <- cmdOpt$selectedGraphics
	
	userOptions$isDispExpDesign <- !regexpr("e",cmdOpt$selectedGraphics) > -1
	userOptions$isFdrPlots <- !regexpr("f",cmdOpt$selectedGraphics) > -1
	userOptions$isIntensityDensityPlots <- !regexpr("i",cmdOpt$selectedGraphics) > -1
	#spectralcountsDensityPlots <- !regexpr("s",cmdOpt$selectedGraphics) > -1
	userOptions$isVolcanoPlots <- !regexpr("v",cmdOpt$selectedGraphics) > -1
	#mvaPlots <- !regexpr("m",cmdOpt$selectedGraphics) > -1
	userOptions$isCorrelationPlots <- !regexpr("c",cmdOpt$selectedGraphics) > -1
	userOptions$isHClustPlot <- !regexpr("h",cmdOpt$selectedGraphics) > -1
	userOptions$isDeFdrPlot <- !regexpr("d",cmdOpt$selectedGraphics) > -1
	#lodsDistributionPlots <- regexpr("l",cmdOpt$selectedGraphics) > -1
	
	# graphics end
	
	### test run
	userOptions$test <- cmdOpt$test
	
	### TEST NON-REQUIRED OPTIONS END
	
	return(userOptions)
}


#print(userOptions)

### SET USER OPTIONS END
