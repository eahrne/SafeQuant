# TODO: Add comment
# 
# Author: erikahrne
###############################################################################
### load / source

##@TEMP

source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Targeted.R")

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

#skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX1b_WB_dilution_series_26112014_Transition Results_plus_background_forErik_filtered_2.csv"
#skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX2b_WB_dilution_series_26112014_Transition Results_plus_background_forErik_final_4.csv"
#skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX3b_WB_dilution_series_26112014_Transition Results_for_Erik_final.csv"

#skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX3b_WB_dilution_series_26112014_Transition Results_for_Erik_final_blank25.csv"
#pdfFile <- "/Users/erikahrne/tmp/NRX3b_WB_dilution_series_26112014_Transition Results_for_Erik_final_blank25_LOW.pdf"

#skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX3b_WB_dilution_series_26112014_Transition Results_for_Erik_final_blank25.csv"
#pdfFile <- "/Users/erikahrne/tmp/test.pdf"

skylineExportFile <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/targetedMS/NRX1a-3_WB_dilution_series_23012015_Transition Results_2_for_Ricky.csv"
pdfFile <- "/Users/erikahrne/tmp/test_low.pdf"
#pdfFile <- "/Users/erikahrne/tmp/test_blank.pdf"

### parse
skylineData <- read.csv(skylineExportFile,sep=",")

names(skylineData)
head(skylineData)
str(skylineData)

### HACK replace concentration 
if(F){
	concentration <- gsub("B14\\-110..\\_NRX1b\\_1:" ,"",as.character(skylineData$Replicate.Name))
	concentration <- gsub("_WB.*","",concentration)
	concentration <- gsub("B14.*","",concentration)
	concentration <- 1/as.numeric(concentration)
	concentration[is.na(concentration)] <- 0
	
	repConc <- data.frame( rev(c(52.23909835,10.4558762,2.116880375,0.472364548,0.149315246,0)),row.names=unique(concentration))
	concentration <- repConc[as.character(concentration),]
}
if(F){
	concentration <- gsub("B14\\-110..\\_NRX..\\_1:" ,"",as.character(skylineData$Replicate.Name))
	concentration <- gsub("_WB.*","",concentration)
	concentration <- gsub("B14.*","",concentration)
	concentration <- 1/as.numeric(concentration)
	concentration[is.na(concentration)] <- 0
	
	repConc <- data.frame( rev(c(36.61394259,6.915043504,1.307507019,0.279922881,0.090591176,0)),row.names=sort(unique(concentration)))
	concentration <- repConc[as.character(concentration),]
}
if(T){
	#concentration <- gsub("B15\\-01....\\_NRX..\\_WB_GFPl:*_*1:" ,"",as.character(skylineData$Replicate.Name))
	concentration <- gsub("B15\\-01..._.*:" ,"",as.character(skylineData$Replicate.Name))
	concentration <- gsub("_WB.*","",concentration)
	concentration <- gsub("B15.*","",concentration)
	concentration <- gsub("_.*","",concentration)
	concentration <- 1/as.numeric(concentration)
	concentration[is.na(concentration)] <- 0
	
	repConc <- data.frame( rev(c(40.86970737,8.135769725	,1.58482937,	0.383470541	,0.161746873,0)),row.names=sort(unique(concentration)))
	concentration <- repConc[as.character(concentration),]
	
}

data.frame(concentration,skylineData$Replicate.Name)

skylineData <- cbind(skylineData,concentration)
skylineData <- skylineData[!grepl("N",skylineData$Area),]
skylineData$Area <- as.numeric(as.character(skylineData$Area))

### rollUp -> create ExpressionSet
nbReplicates <- 3
peptideChargeMzStateConc <-  paste(skylineData$Peptide.Sequence,skylineData$Precursor.Charge,round(skylineData$Precursor.Mz),skylineData$concentration,sep="_")
uniquePeptideChargeStateConc <-  unique(peptideChargeMzStateConc)

expressionMatrix <- data.frame(matrix(nrow=length(unique(unique(peptideChargeMzStateConc))),ncol=nbReplicates),row.names=unique(peptideChargeMzStateConc))
featureAnnotations <- data.frame()
for(pcsc in uniquePeptideChargeStateConc){
	
	peptideChargeMzStateConcSubset <- skylineData[peptideChargeMzStateConc %in% pcsc,]
	
	### add feature data
	c <- 1
	for(rep in unique(peptideChargeMzStateConcSubset$Replicate.Name)){
		expressionMatrix[pcsc,c] <-sum(peptideChargeMzStateConcSubset[peptideChargeMzStateConcSubset$Replicate.Name %in% rep,]$Area)
			c <- c+1
	}

	### add assay data
	featureAnnotations  <- rbind(featureAnnotations,peptideChargeMzStateConcSubset[1,c(1,2,4,5,9,13)])
	
}

featureAnnotations <-cbind(featureAnnotations,dilutionCurveId=paste(featureAnnotations$Peptide.Sequence,featureAnnotations$Precursor.Charge,round(featureAnnotations$Precursor.Mz),sep="_")) 

expressionMatrix[expressionMatrix == 0] <- NA 
row.names(featureAnnotations) <- unique(peptideChargeMzStateConc)	

expDesign <- data.frame(condition=rep(1,3),isControl=rep(T,3))
names(expressionMatrix) <- rownames(expDesign)

### rollUp -> create ExpressionSet END

esetCalibCurve <- createExpressionDataset(expressionMatrix=as.matrix(expressionMatrix),expDesign=expDesign,featureAnnotations=featureAnnotations)

### INIT END



testCreateCalibrationCurve <- function(){
	
	cat(" --- testCreateCalibrationCurve --- \n")
	
	calibCurve <- createCalibrationCurve(esetCalibCurve[fData(esetCalibCurve)$dilutionCurveId == unique(fData(esetCalibCurve)$dilutionCurveId)[2], ])
	
	stopifnot(round(calibCurve$lod,2) == 0.17)
	
	cat(" --- testCreateCalibrationCurve: PASS ALL TEST  --- \n")
	
}

### RUN TESTS

#testCreateCalibrationCurve()

### RUN TESTS END

### GRAPHICS

#pdf("/Users/erikahrne/tmp/NRX1b_WB_dilution_series_26112014_Transition Results_plus_background_forErik_filtered_2_DILUTION_CURVES_blank.pdf")
#pdf("/Users/erikahrne/tmp/NRX1b_WB_dilution_series_26112014_Transition Results_plus_background_forErik_filtered_2_DILUTION_CURVES_lowconc.pdf")

pdf(pdfFile)

par(mfrow=c(1,1))

for(peptide in unique(fData(esetCalibCurve)$Peptide.Sequence)){
		
	### HACK plot heavy only
	idx <- sort(unique(fData(esetCalibCurve)$dilutionCurveId[fData(esetCalibCurve)$Peptide.Sequence %in% peptide]))[2]
	#print(length(sort(unique(fData(esetCalibCurve)$dilutionCurveId[fData(esetCalibCurve)$Peptide.Sequence %in% peptide]))))
	
	#plot(createCalibrationCurve(esetCalibCurve[as.character(fData(esetCalibCurve)$dilutionCurveId) == idx, ],method="blank"), xlab="Concentration (fmol/ul)")
	plot(createCalibrationCurve(esetCalibCurve[as.character(fData(esetCalibCurve)$dilutionCurveId) == idx, ],method="low"), xlab="Concentration (fmol/ul)")
}

graphics.off()


#esetCalib <- esetCalibCurve[as.character(fData(esetCalibCurve)$dilutionCurveId) == as.character(fData(esetCalibCurve)$dilutionCurveId)[75], ]
#calibCurve <- createCalibrationCurve(d,method="low")
#
#exprs(esetCalib)



