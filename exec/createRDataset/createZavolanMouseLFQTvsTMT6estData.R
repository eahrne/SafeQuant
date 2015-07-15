# TODO: Add comment
# 
# Author: erikahrne
###############################################################################
#@TEMP
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")

### PARAMS

## MOUSE
ms1File <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Quant_Ratios_Test_191213/For_Erik_Final/peptides_mouse_FILTERED.csv"
tmtFile <- '/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Quant_Ratios_Test_191213/For_Erik_Final/Raw Data Report for Q+_Zavolan-GM-Mouse-Hepatocytes-TMT-OGE-Final-020913_0%_isotop_Corr_050214.xls'
expDesignTMTSixPlex <- data.frame(condition=paste("cond",sort(rep(c(1,2),2)),sep="_"),isControl=sort(rep(c(T,F),2),decreasing=T) )
pdReportFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/LFQ_vs_TMT/Feature Extractor.txt"

peptidesPerProt <- 2

fdrCutOff <- 0.01
pMassErrorTol <- c(-3,8)

### PARAMS END

### load LFQ MS1 data
esetLFQSpectrum <- parseProgenesisPeptideCsv(file=ms1File,expDesign=getExpDesignProgenesisCsv(ms1File))
# no missing values allowed
esetLFQSpectrum <- esetLFQSpectrum[as.vector(!is.na(apply(exprs(esetLFQSpectrum),1,sum))),]
esetLFQSpectrum <- esetLFQSpectrum[(fData(esetLFQSpectrum)$pMassError >=  pMassErrorTol[1] )  &  (fData(esetLFQSpectrum)$pMassError <=  pMassErrorTol[2]) ,]
esetLFQSpectrum <- addIdQvalues(esetLFQSpectrum)
esetLFQSpectrum <- esetLFQSpectrum[fData(esetLFQSpectrum)$idQValue <=  fdrCutOff,]

### load TMT data
esetTMTSpectrum <- parseScaffoldRawFile(tmtFile,expDesign=expDesignTMTSixPlex, isPurityCorrect=T)
esetTMTSpectrum <- esetTMTSpectrum[!isDecoy(fData(esetTMTSpectrum)$proteinName),]
esetTMTSpectrum <- esetTMTSpectrum[!isCon(fData(esetTMTSpectrum)$proteinName),]
# no missing values allowed
esetTMTSpectrum <- esetTMTSpectrum[as.vector(!is.na(apply(exprs(esetTMTSpectrum),1,sum))),]
sharedPeptides <- intersect(fData(esetLFQSpectrum)$peptide,fData(esetTMTSpectrum)$peptide)

### keep only entries with shared (on peptide level) between MS1 and TMT
esetLFQSpectrum <- esetLFQSpectrum[fData(esetLFQSpectrum)$peptide %in% sharedPeptides,]
esetTMTSpectrum <- esetTMTSpectrum[fData(esetTMTSpectrum)$peptide %in% sharedPeptides,]

#### ROLL UP

### TMT

### parse proteome discoverer report file
# scanNb format "A14-08007.23464X"
scanNb <- paste(gsub("\\.[0-9]{1,6}\\.[0-9]{1,2}$","",fData(esetTMTSpectrum)$spectrumName),"X",sep="")
pdReport <- read.csv(pdReportFile,sep="\t")
rownames(pdReport) <- paste(gsub("\\.raw","",pdReport$Spectrum.File),".",pdReport$First.Scan,"X",sep="") 

### add selected columns to fData
pdAddedColumns <- data.frame(pdReport[scanNb,c("Precursor.Intensity","Isolation.Interference....","Ion.Inject.Time..ms.")])
names(pdAddedColumns) <- c("ms1Int","interference","injectionTime")
fData(esetTMTSpectrum) <- cbind(fData(esetTMTSpectrum),pdAddedColumns)

### calculate mz frequency
mzFreqTable <- table(round(fData(esetTMTSpectrum)$mz)) / sum(table(round(fData(esetTMTSpectrum)$mz)))
mzFreq 	<-  mzFreqTable[as.character(round(fData(esetTMTSpectrum)$mz))]
fData(esetTMTSpectrum) <- cbind(fData(esetTMTSpectrum),mzFreq)

### replace 0's with NA
fData(esetTMTSpectrum)$ms1Int[fData(esetTMTSpectrum)$ms1Int == 0] <- NA
fData(esetTMTSpectrum)$interference[fData(esetTMTSpectrum)$interference == 0] <- NA

esetTMTPeptide <- rollUp(esetTMTSpectrum, featureDataColumnName= c("peptide"), method=c("top1"),isProgressBar=T ) 
esetTMTProtein <- rollUp(esetTMTPeptide, featureDataColumnName= c("proteinName"), method=c("sum"),isProgressBar=T ) 


### LFQ

### same protein (inference) annotations for TMT and LFQ
fData(esetLFQSpectrum)$proteinName <- fData(esetTMTPeptide[as.character(fData(esetLFQSpectrum)$peptide),])$proteinName

esetLFQPeptide <- rollUp(esetLFQSpectrum, featureDataColumnName= c("peptide","ptm"), method=c("top1"),isProgressBar=T ) 
### discard PEPTIDE_ -> PEPTIDE
rownames(esetLFQPeptide) <- gsub("_","",rownames(esetLFQPeptide))

### select for no modified data and peptides with 
esetLFQPeptide <- esetLFQPeptide[nchar(as.character(fData(esetLFQPeptide)$ptm)) == 0,]
esetLFQPeptide <- esetLFQPeptide[!isCon(fData(esetLFQPeptide)$proteinName),]
esetLFQPeptide <- addIdQvalues(esetLFQPeptide)
esetLFQPeptide <- esetLFQPeptide[fData(esetLFQPeptide)$idQValue <=  fdrCutOff,]
esetLFQPeptide <- esetLFQPeptide[!isDecoy(fData(esetLFQPeptide)$proteinName),]

esetLFQProtein <- rollUp(esetLFQPeptide, featureDataColumnName= c("proteinName"), method=c("sum"),isProgressBar=T ) 

### require x peptides per protein
okLFQ <- fData(esetLFQProtein)$nbRolledFeatures >= peptidesPerProt
esetLFQProtein <- esetLFQProtein[okLFQ,]



### FILTER PROTEIN ROLL-UP

okTMT <- fData(esetTMTProtein)$nbRolledFeatures >= peptidesPerProt
esetTMTProtein <- esetTMTProtein[okTMT,]





save(esetLFQSpectrum,esetLFQPeptide,esetLFQProtein,esetTMTSpectrum,esetTMTProtein,esetTMTPeptide,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/zavolanMouseLFQvsTMT6.rda")


print("DONE")
