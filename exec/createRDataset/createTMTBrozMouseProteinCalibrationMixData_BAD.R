# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

#This is 6-plex TMT with the following ratios:
#		
# Phosphorlyase B: 4:1
# Conalbmuin: 1:4
# Ovalbumin: 2:1
# beta-casein: 1:2
# beta-galactosidase 2:1
# alpha-lactablumin: 1:2
#
#Pairs are always 126/127, 128/129 and 130/131. The concentration is highest for 128/129, 5x less in 126/127 and 25x less in 130/131.

################## CONFLICTING PEPTIDES ##################

#sp|P00722|BGAL_ECOLI	Beta-galactosidase OS=Escherichia coli (strain K12) GN=lacZ PE=1 SV=2	59	59	0	14.3898305085	116391.100660001	1024
#sp|P00489|PYGM_RABIT	Glycogen phosphorylase, muscle form OS=Oryctolagus cuniculus GN=PYGM PE=1 SV=3	65	65	51	10.6	97209.8471700004	843
#sp|P02789|TRFE_CHICK	Ovotransferrin OS=Gallus gallus PE=1 SV=2	50	50	3	12.1	77708.3879200003	705
#sp|P01012|OVAL_CHICK	Ovalbumin OS=Gallus gallus GN=SERPINB14 PE=1 SV=2	22	22	3	15.9545454545	42835.49196	386
#tr|B6V3I5|B6V3I5_BOVIN	Alpha-lactalbumin OS=Bos taurus PE=2 SV=1	7	7	0	12	16217.89797	142
#sp|P02666|CASB_BOVIN	Beta-casein OS=Bos taurus GN=CSN2 PE=1 SV=2	9	8	0	11.7777777778	25073.24305	224
#sp|P02662|CASA1_BOVIN	Alpha-S1-casein OS=Bos taurus GN=CSN1S1 PE=1 SV=2	13	12	0	11.9230769231	24495.41514	214
#sp|P02663|CASA2_BOVIN	Alpha-S2-casein OS=Bos taurus GN=CSN1S2 PE=1 SV=2	17	17	2	11.1764705882	25984.25585	222
#
#
# RABBIT
# sp|P00489|PYGM_RABIT	Glycogen phosphorylase, muscle form OS=Oryctolagus cuniculus GN=PYGM PE=1 SV=3	65	65	51	10.6	97209.8471700004	843
#MSRPLSDQEKRKQISVRGLAGVENVTELKKNFNRHLHFTLVKDRNVATPRDYYFALAHTVRDHLVGRWIRTQQHYYEKDPKRIYYLSLEFYMGRTLQNTMVNLALENACDEATYQLGLDMEELEEIEEDAGLGNGGLGRLAACFLDSMATLGLAAYGYGIRYEFGIFNQKICGGWQMEEADDWLRYGNPWEKARPEFTLPVHFYGRVEHTSQGAKWVDTQVVLAMPYDTPVPGYRNNVVNTMRLWSAKAPNDFNLKDFNVGGYIQAVLDRNLAENISRVLYPNDNFFEGKELRLKQEYFVVAATLQDIIRRFKSSKFGCRDPVRTNFDAFPDKVAIQLNDTHPSLAIPELMRVLVDLERLDWDKAWEVTVKTCAYTNHTVLPEALERWPVHLLETLLPRHLQIIYEINQRFLNRVAAAFPGDVDRLRRMSLVEEGAVKRINMAHLCIAGSHAVNGVARIHSEILKKTIFKDFYELEPHKFQNKTNGITPRRWLVLCNPGLAEIIAERIGEEYISDLDQLRKLLSYVDDEAFIRDVAKVKQENKLKFAAYLEREYKVHINPNSLFDVQVKRIHEYKRQLLNCLHVITLYNRIKKEPNKFVVPRTVMIGGKAAPGYHMAKMIIKLITAIGDVVNHDPVVGDRLRVIFLENYRVSLAEKVIPAADLSEQISTAGTEASGTGNMKFMLNGALTIGTMDGANVEMAEEAGEENFFIFGMRVEDVDRLDQRGYNAQEYYDRIPELRQIIEQLSSGFFSPKQPDLFKDIVNMLMHHDRFKVFADYEEYVKCQERVSALYKNPREWTRMVIRNIATSGKFSSDRTIAQYAREIWGVEPSRQRLPAPDEKIP
#
# MOUSE
# sp|Q9WUB3|PYGM_MOUSE	Glycogen phosphorylase, muscle form OS=Mus musculus GN=Pygm PE=1 SV=3	65	65	51	10.6	97206.7053600004	842
#MSRPLSDQDKRKQISVRGLAGVENVSELKKNFNRHLHFTLVKDRNVATPRDYYFALAHTVRDHLVGRWIRTQQHYYEKDPKRIYYLSLEFYMGRTLQNTMVNLALENACDEATYQLGLDMEELEEIEEDAGLGNGGLGRLAACFLDSMATLGLAAYGYGIRYEFGIFNQKICGGWQMEEADDWLRYGNPWEKARPEFTLPVHFYGRVEHTSQGAKWVDTQVVLAMPYDTPVPGYRNNVVNTMRLWSAKAPNDFNLKDFNVGGYIQAVLDRNLAENISRVLYPNDNFFEGKELRLKQEYFVVAATLQDIIRRFKSSKFGSRDPVRTNFDAFPDKVAIQLNDTHPSLAIPELMRILVDLERLDWDKAWDVTVKTCAYTNHTVLPEALERWPVHLMETLLPRHLQIIYEINQRFLNRVAAAFPGDVDRLRRMSLVEEGAVKRINMAHLCIAGSHAVNGVARIHSEILKKTIFKDFYELEPHKFQNKTNGITPRRWLVLCNPGLAEVIAERIGEDYISDLDQLRKLLSYVDDEAFIRDVAKVKQENKLKFSAYLEREYKVHINPNSLFDVQVKRIHEYKRQLLNCLHIITLYNRIKREPNRFMVPRTIMIGGKAAPGYHMAKMIIKLITAIGDVVNHDPAVGDRLRVIFLENYRVSLAEKVIPAADLSEQISTAGTEASGTGNMKFMLNGALTIGTMDGANVEMAEEAGEENFFIFGMRVEDVERLDQRGYNAQEYYDRIPELRQIIEQLSSGFFSPKQPDLFKDIVNMLMHHDRFKVFADYEEYIKCQDKVSELYKNPREWTRMVIRNIATSGKFSSDRTIAQYAREIWGVEPSRQRLPAPDEKI
#


source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")



scaffoldRawDataFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Protmix_111114/Scaffold/Raw Data Report for roland_protein_calMix.xls"
#pdReportFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Protmix_111114/"
conflictingPeptideFile <- "/Volumes/pcf01$/Schmidt_Group/Alex/LC-MS-Tests/TMT_Protmix_111114/Database/s_mouse_tmt_d_DIGEST_conflictingPeptides.csv"

### READ DATA
### parse scaffold file

### parse scaffold file
expDesignTMTSixPlex <- data.frame(condition=paste("cond",c(1:6),sep="_"),isControl=rep(F,6) )
expDesignTMTSixPlex$isControl[c(1)] <- T

esetTMT6Spectrum <- parseScaffoldRawFile(scaffoldRawDataFile, expDesign=expDesignTMTSixPlex, isPurityCorrect=F)

#filter
esetTMT6Spectrum <- esetTMT6Spectrum[!isDecoy(fData(esetTMT6Spectrum)$proteinName),]  
esetTMT6Spectrum <- esetTMT6Spectrum[!isCon(fData(esetTMT6Spectrum)$proteinName),]  
esetTMT6Spectrum <- esetTMT6Spectrum[!(fData(esetTMT6Spectrum)$peptide %in% as.character(read.csv(conflictingPeptideFile)[,1])),]

### filter based on number of peptides per protein
minPepPerProt <- 2
keepACs <- names(table(fData(esetTMT6Spectrum)$proteinName)[table(fData(esetTMT6Spectrum)$proteinName) >= minPepPerProt ])
esetTMT6Spectrum <- esetTMT6Spectrum[(fData(esetTMT6Spectrum)$proteinName %in% keepACs),] 

### parse proteome discoverer report file
# scanNb format "A14-08007.23464X"
if(F){
	scanNb <- paste(gsub("\\.[0-9]{1,6}\\.[0-9]{1,2}$","",fData(esetTMT6Spectrum)$spectrumName),"X",sep="")
	pdReport <- read.csv(pdReportFile,sep="\t")
	rownames(pdReport) <- paste(gsub("\\.raw","",pdReport$Spectrum.File),".",pdReport$First.Scan,"X",sep="") 
	rownames(pdReport) <- gsub("_Ratio","",rownames(pdReport)) 
	
	### add selected columns to fData
	pdAddedColumns <- data.frame(pdReport[scanNb,c("Precursor.Intensity","Isolation.Interference....","Ion.Inject.Time..ms.","Precursor.Area")])
	names(pdAddedColumns) <- c("ms1Int","interference","injectionTime","ms1Area")
	fData(esetTMT6Spectrum) <- cbind(fData(esetTMT6Spectrum),pdAddedColumns)
}


#esetTMT6Peptide <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("peptide"),method="sum",isProgressBar=T) 
#fData(esetTMT6Peptide)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Peptide)$proteinName ) 

esetTMT6Protein <- rollUp(esetTMT6Spectrum,featureDataColumnName= c("proteinName"),method="sum",isProgressBar=T) 
fData(esetTMT6Protein)$isNormAnchor <-  grepl("HUMAN",fData(esetTMT6Protein)$proteinName ) 


#save(esetTMT6Spectrum,esetTMT6Peptide,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )
save(esetTMT6Spectrum,esetTMT6Protein,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/brozProeinCalMixTMT6.rda" )
#save(esetTMT6Spectrum,file="/Users/erikahrne/dev/R/workspace/SafeQuant/data/proteomeMixTMT6.rda" )

print("DONE")
