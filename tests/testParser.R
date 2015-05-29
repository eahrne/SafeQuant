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
	
	cat(" --- testParseProgenesisPeptideMeasurementCsv --- \n")
	eset <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementFile1))
	stopifnot(ncol(exprs(eset)) == 1)
	stopifnot(nrow(exprs(eset)) == 9840)
	
	stopifnot(sum(grepl(";",fData(eset)$proteinName)) == 2)
	stopifnot("sp|Q9Y6H1|CHCH2_HUMAN;sp|Q5T1J5|CHCH9_HUMAN" %in% fData(eset)$proteinName)
	
	cat(" --- testParseProgenesisPeptideMeasurementCsv: PASS ALL TEST  --- \n")
	
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
### TESTS END

