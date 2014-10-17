#!/usr/bin/Rscript

# TODO: Run Safequant
# 
# Author: erikahrne
###############################################################################


# TEST FILE
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/peptides1_FILTERED.csv /Volumes/pcf01\$/Schmidt_Group/Databases/SwissProt_Databases/s_human_d_201405.fasta
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/proteins1.csv
# /Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/TMT_6-Plex_Scaffold_Raw_Export_Example.xls

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

### PARSERS

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
		cat("INFO: PARSING PROGENESIS PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisProteinCsv(file=userOptions$inputFile,expDesign=expDesign)
		
	}else{ #"ProgenesisPeptide"
		cat("INFO: PARSING PROGENESIS FEATURE EXPORT FILE ",userOptions$inputFile, "\n" )
		sqaMethod <- c("rt","naRep")
		eset <- parseProgenesisPeptideCsv(file=userOptions$inputFile,expDesign=expDesign)
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

# parse .fasta file
if(!is.na(userOptions$proteinFastaFile)){
	cat("INFO: PARSING PROTEIN SEQUENCE DB ",userOptions$proteinFastaFile, "\n" )
	### read protein db
	proteinDB <- read.fasta(userOptions$proteinFastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
}

### PARSERS END

if(userOptions$verbose)  print(eset)
if(userOptions$verbose)  print(pData(eset))
if(userOptions$verbose) print(names(fData(eset)))

#### CREATE FEATURE DATA AND FILTER (pre-rollup)

# generic
filter <- data.frame(
		isDecoy(fData(eset)$proteinName) # decoy
		,isCon(fData(eset)$proteinName)	# contaminants	
		,!(grepl(userOptions$selectedProteinName,fData(eset)$proteinName,ignore.case=T)) # protein ac
)

if("idScore" %in% names(fData(eset))){
### applicable to Progenesis Exports	
	
	# add id-level qValues
	eset <- addIdQvalues(eset)
	filter <- cbind(filter,fData(eset)$idQValue > userOptions$fdrCutoff)
}	

if("pMassError" %in% names(fData(eset))){
### applicable to Progenesis feature Exports	

	filter <- cbind(filter, 
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
	filter <- cbind(filter,!(grepl(userOptions$selectedModifName,fData(eset)$ptm,ignore.case=T)))
}

if("nbPeptides" %in% names(fData(eset))){
### applicable to Progenesis Protein Exports (or other upon Protein level roll-up)
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

# no roll-up
#	- progenesis protein

# roll-up peptide-ptm level
#	- progenesis peptide

# roll-up top3
#	- progenesis peptide

# roll-up protein level
#	- scaffold  	
#	- progenesis peptide 

#getIBAQEset
#	- require protein level eset	


### EXPRESSION ANALYSIS

### GRAPHICS

pdf(userOptions$pdfFile)

.qcPlots(sqa$eset,selection=1:6 )
.qcPlots(eset,selection=7 )

cat("INFO: CREATED FILE ", userOptions$pdfFile,"\n")

graphics.off()

### GRAPHICS END

#print(userOptions)

### TSV EXPORT

write.table(cbind(exprs(eset),fData(sqa$eset)),file=userOptions$tsvFilePath,sep="\t", row.names=F)

cat("INFO: CREATED FILE ", userOptions$tsvFilePath,"\n")	
	
### TSV EXPORT END

#print(sqa$qValue)



