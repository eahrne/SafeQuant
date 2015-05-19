# TODO: Add comment
# 
# Author: erikahrne
###############################################################################



source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/initTestSession.R")

### TEST FUNCTIONS

################## TMT ##################

testGetSkipLineNb <- function(){
	cat(" --- testGetSkipLineNb --- \n")
	stopifnot(40 ==  .getSkipLineNb(tmt6PlexRawTestFile))
	stopifnot(40 ==  .getSkipLineNb(tmt10PlexRawTestFile))
	cat(" --- testGetSkipLineNb: PASS ALL TEST --- \n")	
}

testGetNbPlex <- function(){
	cat(" --- testGetNbPlex --- \n")
	stopifnot(10 ==.getNbPlex(tmt10PlexRawTestFile))
	stopifnot(6 ==.getNbPlex(tmt6PlexRawTestFile))
	cat(" --- testGetNbPlex: PASS ALL TEST --- \n")	
}



################## LFQ ##################

testGetProgenesisCsvExpressionColIndices <- function(){
	
	cat(" --- testGetProgenesisCsvProteinIntColIndices --- \n")
	
	stopifnot(sum(14:17 %in% .getProgenesisCsvExpressionColIndices(progenesisProteinCsvFile1) ) == 4 )
	
	stopifnot(sum(31:48 %in% .getProgenesisCsvExpressionColIndices(progenesisFeatureCsvFile1) ) == 18 )
	
	cat(" --- testGetProgenesisCsvProteinIntColIndices: PASS ALL TEST --- \n")	
	
}

testGetExpDesignProgenesisCsv <- function(){
	
	cat(" --- testGetExpDesignProgenesisCsv --- \n")
	
	stopifnot(4 == nrow( getExpDesignProgenesisCsv(progenesisProteinCsvFile1) )) 
	
	stopifnot(18 == nrow(getExpDesignProgenesisCsv(progenesisFeatureCsvFile1) ))
	
	cat(" --- testGetExpDesignProgenesisCsv: PASS ALL TEST  --- \n")
	
}


testParseProgenesisProteinCsv <- function(){
	
	cat(" --- testParseProgenesisProteinCsv --- \n")
	
	eset <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisProteinCsvFile1))
	stopifnot(1917 == nrow(eset))
	stopifnot(6 == ncol(fData(eset)))
	stopifnot(4 == nrow(pData(eset)))
	
	expDesign <- getExpDesignProgenesisCsv(progenesisProteinCsvFile1)[c(3,2),]
	eset2 <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=expDesign)
	
	stopifnot(colnames(eset2) == row.names(pData(eset))[c(3,2)]   )
	### change expDesign
	
	cat(" --- testParseProgenesisProteinCsv: PASS ALL TEST  --- \n")
}

testParseProgenesisFeatureCsv <- function(){
	
	cat(" --- testParseProgenesisFeatureCsv --- \n")
	
	eset <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	stopifnot(496 == nrow(eset))
	names(fData(eset))
	stopifnot(13 == ncol(fData(eset)))
	stopifnot(18 == nrow(pData(eset)))
	
#	expDesign <- getExpDesignProgenesisCsv(progenesisFeatureCsvFile1)[c(3,4:8),]
#	eset2 <- parseProgenesisPeptideCsv(file=progenesisFeatureCsvFile1,expDesign=expDesign)
#	
#	pData(eset2)
#	pData(eset)
	
	cat(" --- testParseProgenesisFeatureCsv: PASS ALL TEST  --- \n")
}

testParseScaffoldRawFile <- function(){
	
	
	cat(" --- testParseScaffoldRawFile --- \n")
	
	# 6 plex	
	
	expDesignTMTSixPlex <- data.frame(condition=sort(rep(c(1,2),3)),isControl=sort(rep(c(T,F),3),decreasing=T) )
	eset <- parseScaffoldRawFile(tmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex)
	stopifnot(329 ==  nrow(eset))
	stopifnot("sp|P38711|RS27B_YEAST,sp|P35997|RS27_YEAST" ==  fData(eset)$proteinName[1])
	stopifnot("sp|P38711|RS27B_YEAST" == fData(parseScaffoldRawFile(tmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex,keepFirstAcOnly=T))$proteinName[1])
	
	expDesignTMTSixPlex <- data.frame(condition=sort(rep(c(1,2),3)),isControl=sort(rep(c(T,F),3),decreasing=T) )
	rownames(expDesignTMTSixPlex) <- rev(1:6)
	eset2 <- parseScaffoldRawFile(tmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex[1:5,])
	stopifnot(exprs(eset)[1,6] == exprs(eset2)[1,1])
	
	# 10 plex	
	expDesignTMTTenPlex <- data.frame(condition=sort(rep(c(1:5),2)),isControl=c(T,T,rep(F,8)) )
	eset3 <- parseScaffoldRawFile(tmt10PlexRawTestFile,expDesign=expDesignTMTTenPlex, keepFirstAcOnly=T)
	stopifnot(8875 ==  nrow(eset3))
	cat(" --- testParseScaffoldRawFile: PASS ALL TEST --- \n")
	
}

testGetFileType <- function(){
	
	cat(" --- testGetFileType --- \n")
	
	stopifnot(.getFileType(tmt6PlexRawTestFile) == "ScaffoldTMT")
	stopifnot(.getFileType(tmt10PlexRawTestFile) == "ScaffoldTMT")
	stopifnot(.getFileType(progenesisFeatureCsvFile1) == "ProgenesisFeature")
	stopifnot(.getFileType(progenesisProteinCsvFile1) == "ProgenesisProtein")
	stopifnot(.getFileType(progenesisPeptideMeasurementFile1) == "ProgenesisPeptide")
	stopifnot(.getFileType(maxQuantProteinFileTxt) == "MaxQuantProteinGroup")

	cat(" --- testGetFileType: PASS ALL TEST --- \n")
	
}

testParseProgenesisPeptideMeasurementCsv <- function(){
	
	### cannot use "peptide measurement exports" as feature are mapped to peptides. The same feature is listed multiple times
	
	cat(" --- testParseProgenesisPeptideMeasurementCsv --- \n")
	eset <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementFile1))
	stopifnot(ncol(exprs(eset)) == 1)
	stopifnot(nrow(exprs(eset)) == 13035)
	cat(" --- testParseProgenesisPeptideMeasurementCsv: PASS ALL TEST  --- \n")
	
	
	peptideMultipleUnsorted <- "AAIDWFDGK"
	fData(eset)[fData(eset)$peptide == 	peptideMultipleUnsorted,]
	
	which(fData(eset)$peptide == 	peptideMultipleUnsorted)
	
	### sort protein accessions
#	fData(eset)$proteinName <- as.character(fData(eset)$proteinName)
#	fData(eset)$proteinName <- as.vector(unlist(lapply(fData(eset)$proteinName, function(prot){
#	
#				prot <- paste(sort(as.vector(unlist(strsplit(prot,";")))),collapse=";")
#				return(prot)
#				
#			}  )))
	
	
}


testParseMaxQuantProteinGroupTxt <- function(){
	
	cat(" --- testParseMaxQuantProteinGroupTxt:  --- \n")
	
	expDesign <- data.frame(row.names=1:20,isControl=c(rep(T,5),rep(F,15)),condition=sort(rep(paste("condition",1:4,sep="_"),5)))
	eset <- parseMaxQuantProteinGroupTxt(maxQuantProteinFileTxt,expDesign=expDesign, method="auc")
	
	stopifnot(ncol(fData(eset)) == 7)
	stopifnot(ncol(eset) == 20)
	stopifnot( nrow(exprs(parseMaxQuantProteinGroupTxt(maxQuantProteinFileTxt,expDesign=expDesign, method="spc"))) == 3426 )

	cat(" --- testParseMaxQuantProteinGroupTxt: PASS ALL TEST  --- \n")
	
}


### TEST FUNCTIONS END


### TESTS
if(F){
	testGetSkipLineNb()
	testParseScaffoldRawFile()
	testGetNbPlex()
	testGetProgenesisCsvExpressionColIndices()
	testGetExpDesignProgenesisCsv()
	testParseProgenesisProteinCsv()
	testParseProgenesisFeatureCsv()
	testParseProgenesisPeptideMeasurementCsv()
	testParseMaxQuantProteinGroupTxt()
	testGetFileType()
}
### TESTS END

# OCCAMS RAZOR IMPLEMENTATION (HIGHEST SCORING PROTEIN TAKES IT ALL)
# NOT COMPATIBLE WITH STANDARD iBAQ AND TOP3 ABSOLUT QUANTIFICATION! # EXAMPLE: Peptide Measurement Output
# Accession	sp|O43707|ACTN4_HUMAN
# All accessions (for this sequence)	sp|O43707|ACTN4_HUMAN;sp|P12814|ACTN1_HUMAN;sp|P35609|ACTN2_HUMAN;sp|Q08043|ACTN3_HUMAN
# Grouped accessions (for this sequence)	sp|P35609|ACTN2_HUMAN;sp|Q08043|ACTN3_HUMAN
# Shared accessions (for this sequence)	sp|P12814|ACTN1_HUMAN
			
# PROTEIN GROUP: TWO PROTEINS SHARING A SET OF PEPTIDES
# EXAMPLE - A GROUP
# Protein X: PEP-A,PEP-B
# Protein Y: PEP-A,PEP-B

# EXAMPLE: SHARED SUBSET - ONE PROTEIN HAS SPECIFIC PEPTIDE(S)
# Protein X: PEP-A,PEP-B, PEP-C
# Protein Y: PEP-A,PEP-B
# Here all petides should be assigned to Protein X and Protein Y listed as "SUBSET" Protein

# EXAMPLE - SHARED SUBSET - BOTH PROTEINS HAVE SPECIFIC PEPTIDE(S)
# Protein X: PEP-A,PEP-B, PEP-C, PEP-D
# Protein Y: PEP-A,PEP-B, PEP-E
# If Score(X) > score(Y)
# -> PEP-A,PEP-B, PEP-C, PEP-D Assigned to Protein X  
# -> PEP-E Assigned to Protein Y


# EXAMPLE - SHARED SUBSET - AT LEAST TWO PROTEINS HAVE SPECIFIC PEPTIDE(S)
# Protein X: PEP-A,PEP-B, PEP-C, PEP-D
# Protein Y: PEP-A,PEP-B, PEP-E
# Protein Z: PEP-A,PEP-B, PEP-C, PEP-D
# GROUP Protein X and Protein Z
# -> PEP-A,PEP-B, PEP-C, PEP-D Assigned to GROUP Protein X;  Protein Z
# -> PEP-E Assigned to Protein Y


# EXAMPLE - SHARED SUBSET - AT LEAST TWO PROTEINS HAVE SPECIFIC PEPTIDE(S)
# Protein X: PEP-A,PEP-B, PEP-C, PEP-D
# Protein Y: PEP-A,PEP-B, PEP-E
# Protein Z: PEP-A,PEP-B, PEP-E, PEP-F 
# GROUP Protein X and Protein Z
# If ((Score(X) > score(Y)) &  (Score(X) > score(Z)))
# -> PEP-A,PEP-B, PEP-C, PEP-D Assigned to Protein X  
# -> PEP-E, PEP-F Assigned to Protein Z


### EXAMPLE PROTEIN LISTED AS GROUPED AND NOT GROUPE
#SEQUENCE 	VILLGDGGVGK	
#Accession sp|P51151|RAB9A_HUMAN
#Grouped accessions sp|Q9NP90|RAB9B_HUMAN
#
#SEQUENCE	DATNVAAAFEEAVR	
#Accession sp|P51151|RAB9A_HUMAN
#Grouped accessions 
#
# Here all peptides mapped to sp|Q9NP90|RAB9B_HUMAN also map to sp|P51151|RAB9A_HUMAN, while additional peptides are mapped to sp|P51151|RAB9A_HUMAN 

library(dplyr)
library(data.table)

#file <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements1.csv"
#file <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements2.csv"
file <- "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements3.csv"

### 

#	Feature	Score	Sequence	Accession				All accessions (for this sequence)	Grouped accessions (for this sequence)	Shared accessions (for this sequence)
#	1770	62.91	AALSALESFLK	sp|P78527|PRKDC_HUMAN	sp|P78527|PRKDC_HUMAN		

# A "Grouped accessions (for this sequence)" never appears in the "Accession" column as ProteinScore(GroupedProtein) <= ProteinScore(Accession)
# allGrouped <-  as.vector(unlist(lapply(resDT$"Grouped accessions (for this sequence)",function(t){strsplit(as.character(t),";")})))
# sum(allGrouped %in% as.character(resDT$Accession))

# OCCAMS RAZOR IMPLEMENTATION (HIGHEST SCORING PROTEIN/PROTEIN GROUP TAKES IT ALL)
# NOT COMPATIBLE WITH STANDARD iBAQ AND TOP3 ABSOLUT QUANTIFICATION! 

res <- read.csv(file,skip=2,allowEscapes=T,check.names=F)



#names(resDT)[14] <- "SharedAccessions"
#resDT$SharedAccessions <- as.character(resDT$SharedAccessions)
#bestProteinPerFeatureIdx <-resDT[, list(  selectedResDTIndex  = resDTIndex[order(proteinScore,decreasing=T)[1]], altAccessions= paste(unique(unlist(strsplit(SharedAccessions, ";"))),collapse=";") ), by = key(resDT)]
#bestProteinPerFeatureIdx <-resDT[, list( tmpAC =  Accession[order(proteinScore,decreasing=T)[1]],selectedResDTIndex  = resDTIndex[order(proteinScore,decreasing=T)[1]], altAccessions= paste(unique(unlist(strsplit(SharedAccessions, ";"))),collapse=";") ), by = key(resDT)]



# DONE

#resDT[resDT$Feature == 2,]
#featureDT[featureDT$Feature == 2,]
#
#rownames(bestProteinPerFeatureIdx) <- bestProteinPerFeatureIdx$selectedResDTIndex
#featureDT <- cbind(resDT[bestProteinPerFeatureIdx$selectedResDTIndex,],altAccessions=bestProteinPerFeatureIdx[bestProteinPerFeatureIdx$selectedResDTIndex,]$altAccessions,tmpAC= bestProteinPerFeatureIdx[bestProteinPerFeatureIdx$selectedResDTIndex,]$tmpAC)
#data.frame(featureDT$tmpAC,featureDT$Accession)
#
#featureDT[featureDT$Accession %in% leadingACTrueGroup,]
#resDT[resDT$Accession %in% leadingACTrueGroup,]
#
#
#names(bestProteinPerFeatureIdx)
#
#bestProteinPerFeatureIdx$altAccessions[10852]






method="auc"
#file = "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements1.csv"
#file = "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements2.csv"
#file = "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements3.csv"
file = "/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements4.csv"
expDesign= getExpDesignProgenesisCsv(file)






#eset <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementFile1))
#stopifnot(ncol(exprs(eset)) == 1)
#stopifnot(nrow(exprs(eset)) == 13035)

print("DONE")
