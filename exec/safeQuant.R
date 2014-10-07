#!/usr/bin/Rscript

# TODO: Run Safequant
# 
# Author: erikahrne
###############################################################################


# TEST FILE
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/PeptidesSQAnalysis/peptides1_FILTERED.csv
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/ProteinsSQAnalysis/proteins1.csv
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/TMT/TMT_Scaffold_Raw_Export_Example.xls

if(F){
	### 
	library("SafeQuant")
}else{
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
}

VERSION <- 2.01

### USER CMD LINE OPTIONS
userOptions <- getUserOptions(version=VERSION)
### USER CMD LINE OPTIONS END

### SUPRESS WARNINGS
if(!userOptions$verbose){
	options(warn=-1)
}

#### DEPENDENCIES

if(userOptions$verbose) cat("LOADING DEPENDENCIES \n")	

suppressPackageStartupMessages(library("affy", quiet=T))
suppressPackageStartupMessages(library("limma", quiet=T))
suppressPackageStartupMessages(library(gplots, quiet=T)) # volcano plot
suppressPackageStartupMessages(library(seqinr, quiet=T))

### PARSE INPUT FILE

if(userOptions$verbose) cat("PARSING INPUT FILE \n")	

# get file type
fileType <- .getFileType(userOptions$inputFile)

sqaMethod <- c("global","naRep")

### Progenesis Export
if(fileType %in% c("ProgenesisProtein","ProgenesisPeptide")){
	
	# default
	expDesign <- getExpDesignProgenesisCsv(userOptions$inputFile)
	
	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)
	}
	
	if(fileType == "ProgenesisProtein"){
		eset <- parseProgenesisProteinCsv(file=userOptions$inputFile,expDesign=expDesign)
		
		# fdr filter
		# ac filter
		# peptides per protein filter
		# replace missing values
		# normalize
		# stat test
		# (abs. est)
		
	}else{ #"ProgenesisPeptide"
		sqaMethod <- c("rt","naRep")
		eset <- parseProgenesisPeptideCsv(file=userOptions$inputFile,expDesign=expDesign)
		
		# fdr filter
		# ac filter
		# ptm filter
		# (peptides per protein filter)
		# mass accuracy filter
		# replace missing values
		# normalize
		# (roll up)
		# (protein inference)
		# stat test
		# (abs. est)
	}
	
# Scaffold Export (TMT data)
}else if(fileType == "ScaffoldTMT"){
	
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
### PARSE INPUT FILE END

if(userOptions$verbose)  print(eset)
if(userOptions$verbose)  print(pData(eset))

# ProgenesisProtein

# fdr filter
# ac filter
# peptides per protein filter
# replace missing values
# normalize 
# stat test
# (abs. est)
# graphics	
# tsv export

# ProgenesisPeptide"

# fdr filter
# ac filter
# ptm filter
# (peptides per protein filter)
# mass accuracy filter
# replace missing values
# normalize
# (roll up, peptide/protein)
# (protein inference)
# stat test
# (abs. est)
# graphics
# extract ptm motifs
# tsv export

# ScaffoldTMT

# ac filter
# replace missing values
# normalize
# peptides per protein filter !!!
# roll up protein
# stat test
# graphics
# tsv export

if(userOptions$verbose) print(names(fData(eset)))

### CREATE FEATURE DATA

if("idScore" %in% names(fData(eset))){
	eset <- addIdQvalues(eset)
}

### CREATE FEATURE DATA END

#### CREATE FEATURE DATA AND FILTER (pre-rollup)

# generic
filter <- data.frame(
		isDecoy(fData(eset)$proteinName) # decoy
		,isCon(fData(eset)$proteinName)	# contaminants	
		,!(grepl(userOptions$selectedProteinName,fData(eset)$proteinName,ignore.case=T)) # protein ac
)

if("idScore" %in% names(fData(eset))){
	# add id-level qValues
	eset <- addIdQvalues(eset)
	filter <- cbind(filter,fData(eset)$idQValue > userOptions$fdrCutoff)
}	


if("pMassError" %in% names(fData(eset))){
	
	filter <- cbind(filter, 
			(fData(eset)$pMassError < userOptions$precursorMassFilter[1])
					| (fData(eset)$pMassError > userOptions$precursorMassFilter[2]) # precursor mass tolerance
	)
}

if("ptm" %in% names(fData(eset))){
	
	# add motif-X and ptm coordinates
	if(!is.na(userOptions$proteinFastaFile)){
		eset <- .addPTMCoord(eset,userOptions$proteinFastaFile,motifLength=4)
	}
	filter <- cbind(filter,!(grepl(userOptions$selectedModifName,fData(eset)$ptm,ignore.case=T)))
}

if("nbPeptides" %in% names(fData(eset))){
	filter <- cbind(filter,fData(eset)$nbPeptides < userOptions$minNbPeptidesPerProt)
}	

# set pre-rollup filters
eset <- .setFilter(eset,filter=filter)

#print(fData(eset)$isFiltered)

#### CREATE FEATURE DATA AND FILTER END

### SET ANCHOR PROTEINS
fData(eset)$isNormAnchor <- grepl(userOptions$normAC,fData(eset)$proteinName)

if(userOptions$verbose){
	cat("\nNB. ANCHOR PROTEINS: ")
	cat(sum(fData(eset)$isNormAnchor))
	cat("\n")
	print(rownames(eset)[fData(eset)$isNormAnchor])
	cat("\n")
}

### SET ANCHOR PROTEINS END

### EXPRESSION ANALYSIS

if(userOptions$verbose) print(names(fData(eset)))
sqa <- safeQuantAnalysis(eset, method=sqaMethod)

### EXPRESSION ANALYSIS

### GRAPHICS

pdf(userOptions$pdfFile)

.qcPlots(eset)

cat("CREATED FILE ", userOptions$pdfFile,"\n")

dev.off()

### GRAPHICS END

#print(sqa$qValue)



 #print(userOptions)