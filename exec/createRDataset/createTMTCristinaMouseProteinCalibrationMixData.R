# TODO: Add comment
# 
# Author: erikahrne
###############################################################################



#Phosphorlyase B: 4:1
#Conalbmuin: 1:4
#Ovalbumin: 2:1
#beta-casein: 1:2
#beta-galactosidase 2:1
#alpha-lactablumin: 1:2
#
#The mix was spiked into all samples (channel 1&2 is pair, 3&4 ). We spiked the protein mix in 3 different concentrations 20% (1&2), 100% (3&4), 4% (5&6), 20% (7&8) and 100 % (9&10).

################## CONFLICTING PEPTIDES IN MOUSE ##################
#sp|P00722|BGAL_ECOLI	Beta-galactosidase OS=Escherichia coli (strain K12) GN=lacZ PE=1 SV=2	59	59	0	14.3898305085	116391.100660001	1024
#sp|P00489|PYGM_RABIT	Glycogen phosphorylase, muscle form OS=Oryctolagus cuniculus GN=PYGM PE=1 SV=3	65	65	51	10.6	97209.8471700004	843
#sp|P02789|TRFE_CHICK	Ovotransferrin OS=Gallus gallus PE=1 SV=2	50	50	3	12.1	77708.3879200003	705
#sp|P01012|OVAL_CHICK	Ovalbumin OS=Gallus gallus GN=SERPINB14 PE=1 SV=2	22	22	3	15.9545454545	42835.49196	386
#tr|B6V3I5|B6V3I5_BOVIN	Alpha-lactalbumin OS=Bos taurus PE=2 SV=1	7	7	0	12	16217.89797	142
#sp|P02666|CASB_BOVIN	Beta-casein OS=Bos taurus GN=CSN2 PE=1 SV=2	9	8	0	11.7777777778	25073.24305	224
#sp|P02662|CASA1_BOVIN	Alpha-S1-casein OS=Bos taurus GN=CSN1S1 PE=1 SV=2	13	12	0	11.9230769231	24495.41514	214
#sp|P02663|CASA2_BOVIN	Alpha-S2-casein OS=Bos taurus GN=CSN1S2 PE=1 SV=2	17	17	2	11.1764705882	25984.25585	222
#tr|B6V3I5|B6V3I5_BOVIN	Alpha-lactalbumin OS=Bos taurus PE=2 SV=1	7	7	0	12

source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")

scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20141125-134856_TMTexp20/Scaffold/Raw Data Report for Nigg-Cristina-TMTexp20-081214-calprot.xls"
pdReportFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/20141125-134856_TMTexp20/ProteomeDiscoverer/Cristina_Calib_mix.csv"
conflictingPeptideFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/ENigg/CristinaVigano_93/Databases/digestedProteins_conflictingPeptides.csv"

### READ DATA
### parse scaffold file

### parse scaffold file
expDesignTMT10Plex <- data.frame(condition=c("cond_1","cond_1","cond_2","cond_2",paste("cond",c(3:8),sep="_")),isControl=rep(F,10) )
expDesignTMT10Plex$isControl[c(1,2)] <- T

esetTMT10Spectrum <- parseScaffoldRawFile(scaffoldRawDataFile, expDesign=expDesignTMT10Plex, isPurityCorrect=F)

#filter
esetTMT10Spectrum <- esetTMT10Spectrum[!isDecoy(fData(esetTMT10Spectrum)$proteinName),]  
esetTMT10Spectrum <- esetTMT10Spectrum[!isCon(fData(esetTMT10Spectrum)$proteinName),]  
esetTMT10Spectrum <- esetTMT10Spectrum[!(fData(esetTMT10Spectrum)$peptide %in% as.character(read.csv(conflictingPeptideFile)[,1])),]

### filter based on number of peptides per protein
minPepPerProt <- 2
keepACs <- names(table(fData(esetTMT10Spectrum)$proteinName)[table(fData(esetTMT10Spectrum)$proteinName) >= minPepPerProt ])
esetTMT10Spectrum <- esetTMT10Spectrum[(fData(esetTMT10Spectrum)$proteinName %in% keepACs),] 

### parse proteome discoverer report file
# scanNb format "A14-08007.23464X"

scanNb <- paste(gsub("\\.[0-9]{1,6}\\.[0-9]{1,2}$","",fData(esetTMT10Spectrum)$spectrumName),"X",sep="")
pdReport <- read.csv(pdReportFile,sep=",")
rownames(pdReport) <- paste(gsub("\\.raw","",pdReport$Spectrum.File),".",pdReport$First.Scan,"X",sep="") 
rownames(pdReport) <- gsub("_Ratio","",rownames(pdReport)) 
	
### add selected columns to fData
pdAddedColumns <- data.frame(pdReport[scanNb,c("Precursor.Intensity","Isolation.Interference....","Ion.Inject.Time..ms.","Precursor.Area")])
names(pdAddedColumns) <- c("ms1Int","interference","injectionTime","ms1Area")
fData(esetTMT10Spectrum) <- cbind(fData(esetTMT10Spectrum),pdAddedColumns)

#esetTMT6Peptide <- rollUp(esetTMT10Spectrum,featureDataColumnName= c("peptide"),method="sum",isProgressBar=T) 
#fData(esetTMT6Peptide)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Peptide)$proteinName ) 

esetTMT10Protein <- rollUp(esetTMT10Spectrum,featureDataColumnName= c("proteinName"),method="sum",isProgressBar=T) 
fData(esetTMT10Protein)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT10Protein)$proteinName ) 


#save(esetTMT10Spectrum,esetTMT6Peptide,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )
save(esetTMT10Spectrum,esetTMT10Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/cristinaProeinCalMixTMT10.rda" )
#save(esetTMT10Spectrum,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/cristinaProeinCalMixTMT10.rda" )

print("DONE")
