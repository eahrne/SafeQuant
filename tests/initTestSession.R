# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### load / source

##@TEMP
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/UserOptions.R")

#install.packages("/Users/erikahrne/dev/R/workspace/SafeQuant/", repos = NULL, type="source")
#library(SafeQuant)

library("affy")
library("limma")
library(gplots) # volcano plot
library(seqinr)


### INIT

	### VARIOUS TEST FILES

#testDir <- dirname(sys.frame(1)$ofile)
#testDir <- gsub("tests\\/tmp","inst/tests/",testDir)
tmt6PlexRawTestFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new//TMT_6-Plex_Scaffold_Raw_Export_Example.xls"
tmt10PlexRawTestFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new//TMT_10-Plex_Scaffold_Raw_Export_Example.xls"
progenesisProteinCsvFile1 <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/proteins1.csv"
progenesisProteinCsvFile2 <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/proteins2.csv"

progenesisFeatureCsvFile1 <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/peptides1_FILTERED.csv"
progenesisFeatureCsvFile2 <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/peptides2.csv"

fastaFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/new/sp_mouse_160512.decoy.fasta"


### INIT END

## CREATE TEST DATA

set.seed(1234)
nbFeatures <- 900

peptide <- paste("pep",1:nbFeatures,sep="")
peptide[1] <- "VALGDGVQLPPGDYSTTPGGTLFSTTPGGTR"
peptide[2] <- "AQAGLTATDENEDDLGLPPSPGDSSYYQDQVDEFHEAR"

proteinName <- sort(rep(paste("prot",1:(nbFeatures/3),sep=""),3))
proteinName[1:200] <- paste("REV_",proteinName[1:200] ,sep="")
proteinName[1] <- "sp|Q60876|4EBP1_MOUSE"
proteinName[2] <- "sp|Q9JI13|SAS10_MOUSE"

idScore <- rep(0,length(proteinName))
idScore[1:200] <-rnorm(200,10,1)
idScore[c(1:2,201:900)] <- rnorm(702,15,1)

ptm <- rep("",900)
ptm[1] <- "[15] Phospho (ST)|[30] Phospho (ST)"
ptm[2] <- "[20] Phospho (ST)"

pMassError <- c(rnorm(200,0,1.5),rnorm(700,0,0.5))

charge <- round(runif(length(ptm),1.5,3.7))

peptideName <- paste(peptide,ptm)

proteinDescription <- sort(rep(paste("protDescription",1:(nbFeatures/3),sep=""),3))
isNormAnchor <- rep(T,nbFeatures)
isFiltered <- rep(F,nbFeatures)

m <- as.matrix( data.frame(rnorm(nbFeatures,1001),rnorm(nbFeatures,1001),rnorm(nbFeatures,1002),rnorm(nbFeatures,1002),rnorm(nbFeatures,1000),rnorm(nbFeatures,1000)) )
rownames(m) <- peptideName
colnames(m) <- c("A_rep_1","A_rep_2","B_rep_1","B_rep_2","C_rep_1","C_rep_2")

### phenoData: stores expDesign
#condition isControl
#A_rep_1         A     FALSE
#A_rep_2         A     FALSE
#B_rep_1         B     FALSE
#B_rep_2         B     FALSE
#C_rep_1         C      TRUE
#C_rep_2         C      TRUE

expDesign <- data.frame(condition=c("A","A","B","B","C","C"),isControl=c(F,F,F,F,T,T),row.names=colnames(m))
#expDesign <- data.frame(condition=c("A","A","B","B","C","C"),row.names=colnames(m))

featureAnnotations <- data.frame(
		 peptide
		, charge		 	
		,proteinName
		,proteinDescription
		,idScore
		,ptm
		,pMassError
		,isNormAnchor
		,isFiltered
		,row.names=peptideName)

eset <- createExpressionDataset(expressionMatrix=m,expDesign=expDesign,featureAnnotations=featureAnnotations)
sqa <- safeQuantAnalysis(eset)

# ABS. QUANT SIM. DATA 
cpc <- rep(2^(1:5),10)
set.seed(1234)
signal <- rnorm(length(cpc),cpc,cpc/10)
df <- data.frame(cpc =  log10(cpc),signal = log10(signal))
fit <- lm(cpc ~ signal, data=df )


data(proteomeMixLFQ,package="SafeQuant")
data(proteomeMixTMT6,package="SafeQuant")
