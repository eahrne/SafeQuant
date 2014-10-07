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

progenesisFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/LFQ_Stats_Ratio_Test_09012014/LTQ_Ratio_Statistics_Test_120000_5_replicates_130114/proteins.csv"
expDesign <- getExpDesignProgenesisCsv(progenesisFile)

expDesign$isControl <- rev(expDesign$isControl)

esetLFQProtein <- parseProgenesisProteinCsv(file=progenesisFile,expDesign=expDesign)
### set human proteins as anchors
fData(esetLFQProtein)$isNormAnchor <-  regexpr("HUMAN",fData(esetLFQProtein)$proteinName ) > -1


esetLFQProtein <- addIdQvalues(esetLFQProtein)



filter <- data.frame(
		
		fData(addIdQvalues(esetLFQProtein))$idQValue >= 0.01 # id score
		,isDecoy(fData(esetLFQProtein)$proteinName)	# decoy
		,isCon(fData(esetLFQProtein)$proteinName)	# contaminants	

)

esetLFQProtein <- .setFilter(esetLFQProtein, filter=filter)

### discard BART HUMAN "contaminants"
contaminantsFile <-  "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Ratio_Stats_Paper_230714/Ratio_Pilot_6plex_250714/Scaffold/Human_Contaminants_290714.txt"

contaminantsAC <- read.csv(contaminantsFile,sep="\t",skip=2)
contaminantsAC <- unique(contaminantsAC$Accession.Number)
contaminantsAC <- gsub(" .*","",contaminantsAC)
contaminantsAC <- contaminantsAC[grepl("HUMAN",contaminantsAC)]

esetLFQProtein <- esetLFQProtein[!(fData(esetLFQProtein)$proteinName %in% contaminantsAC),]

save(esetLFQProtein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixLFQ.rda" )

print("DONE")
