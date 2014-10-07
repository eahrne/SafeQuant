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



scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Raw Data Report for TMT-6-stats-2-Bart-human-310714.xls"
pdfFile <- "/Users/erikahrne/dev/R/workspace/TMTRatioCorrection/out/Figure_X_statTest_TMT_2.pdf"
rDataTmpFile <- "/Users/erikahrne/dev/R/workspace/TMTRatioCorrection/rData/statTest_TMT_2.rData"

contaminantsFile <-  "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Human_Contaminants_290714.txt"

### READ DATA



contaminantsAC <- read.csv(contaminantsFile,sep="\t",skip=2)
contaminantsAC <- unique(contaminantsAC$Accession.Number)
contaminantsAC <- gsub(" .*","",contaminantsAC)
contaminantsAC <- contaminantsAC[grepl("HUMAN",contaminantsAC)]

### parse scaffold file

expDesignTMTSixPlex <- data.frame(condition=rep(paste("cond",c(1,2),sep="_"),3),isControl=rep(F,6) )
expDesignTMTSixPlex$isControl[c(1,3,5)] <- T
#expDesignTMTSixPlex$isControl[3] <- T

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

esetTMT6Peptide <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("peptide"),method="sum",isProgressBar=T) 
fData(esetTMT6Peptide)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Peptide)$proteinName ) 

esetTMT6Protein <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("proteinName"),method="sum",isProgressBar=T) 
fData(esetTMT6Protein)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Protein)$proteinName ) 


save(esetTMT6Spectrum,esetTMT6Peptide,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )

print("DONE")
