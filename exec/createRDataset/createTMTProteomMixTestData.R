# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")



### ??? scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Raw Data Report for TMT-6-stats-2-Bart-human-310714.xls"

scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Raw Data Report for Copy_of_TMT-6-ratio-pilot-Bart-human-290714.xls"
pdReportFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Proteome_Discoverer/A14-07151_Ratio-01.txt"

contaminantsFile <-  "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Human_Contaminants_290714.txt"

### READ DATA

contaminantsAC <- read.csv(contaminantsFile,sep="\t",skip=2)
contaminantsAC <- unique(contaminantsAC$Accession.Number)
contaminantsAC <- gsub(" .*","",contaminantsAC)
contaminantsAC <- contaminantsAC[grepl("HUMAN",contaminantsAC)]

### parse scaffold file

### parse scaffold file
expDesignTMTSixPlex <- data.frame(condition=paste("cond",c(1,2,3,1,4,5),sep="_"),isControl=rep(F,6) )
expDesignTMTSixPlex$isControl[c(1,4)] <- T

esetTMT6Spectrum <- parseScaffoldRawFile(scaffoldRawDataFile, expDesign=expDesignTMTSixPlex, isPurityCorrect=F)

#filter
esetTMT6Spectrum <- esetTMT6Spectrum[!isDecoy(fData(esetTMT6Spectrum)$proteinName),]  
esetTMT6Spectrum <- esetTMT6Spectrum[!isCon(fData(esetTMT6Spectrum)$proteinName),]  
esetTMT6Spectrum <- esetTMT6Spectrum[!(fData(esetTMT6Spectrum)$proteinName %in% contaminantsAC),]  
### filter based on number of peptides per protein
minPepPerProt <- 2
keepACs <- names(table(fData(esetTMT6Spectrum)$proteinName)[table(fData(esetTMT6Spectrum)$proteinName) >= minPepPerProt ])
esetTMT6Spectrum <- esetTMT6Spectrum[(fData(esetTMT6Spectrum)$proteinName %in% keepACs),] 
fData(esetTMT6Spectrum)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Spectrum)$proteinName ) 

### parse proteome discoverer report file
# scanNb format "A14-08007.23464X"
scanNb <- paste(gsub("\\.[0-9]{1,6}\\.[0-9]{1,2}$","",fData(esetTMT6Spectrum)$spectrumName),"X",sep="")
pdReport <- read.csv(pdReportFile,sep="\t")
rownames(pdReport) <- paste(gsub("\\.raw","",pdReport$Spectrum.File),".",pdReport$First.Scan,"X",sep="") 
rownames(pdReport) <- gsub("_Ratio","",rownames(pdReport)) 

### add selected columns to fData
pdAddedColumns <- data.frame(pdReport[scanNb,c("Precursor.Intensity","Isolation.Interference....","Ion.Inject.Time..ms.","Precursor.Area")])
names(pdAddedColumns) <- c("ms1Int","interference","injectionTime","ms1Area")
fData(esetTMT6Spectrum) <- cbind(fData(esetTMT6Spectrum),pdAddedColumns)

#esetTMT6Peptide <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("peptide"),method="sum",isProgressBar=T) 
#fData(esetTMT6Peptide)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Peptide)$proteinName ) 

esetTMT6Protein <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("proteinName"),method="sum",isProgressBar=T) 
fData(esetTMT6Protein)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Protein)$proteinName ) 


#save(esetTMT6Spectrum,esetTMT6Peptide,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )
save(esetTMT6Spectrum,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )
#save(esetTMT6Spectrum,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )

print("DONE")
