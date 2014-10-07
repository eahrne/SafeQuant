
### INIT

source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/initTestSession.R")

### INIT END



### TEST FUNCTIONS

testExport <-function(){
	
	cat("--- testExport: --- \n")
	sqa <- safeQuantAnalysis(eset)
	export(sqa,file="/Users/erikahrne/tmp/tmp.csv")
	cat("--- testExport: PASS ALL TEST --- \n")
	
}

testSafequantAnalysis <- function(){
	
	cat("--- testSafequnatAnalysis: --- \n")
	sqa <- safeQuantAnalysis(eset)
	stopifnot(sum(c("eset", "cv", "ratio", "pValue","qValue","baselineIntensity")  %in%  names(sqa)) == 6)
	stopifnot(nrow(sqa$pValue) == nrow(sqa$eset) )
	
	### filter our first 10 peptides
	esetFiltered <- eset
	fData(esetFiltered)$isFiltered[1:10] <- rep(T,10)
	sqaFiltered <- safeQuantAnalysis(esetFiltered)
	stopifnot(sum(is.na(sqaFiltered$pValue[,1])) == 10  )
	stopifnot(sum(!is.na(sqaFiltered$pValue[,1])) == (nrow(sqaFiltered$pValue)-10))
	
	# NA replace
	esetNaRep <- eset
	exprs(esetNaRep)[1,1] <- NA
	stopifnot(!is.na(exprs(safeQuantAnalysis(esetNaRep, method=c("naRep"))$eset)[1,1]))
	stopifnot(is.na(exprs(safeQuantAnalysis(esetNaRep, method=c(""))$eset)[1,1]))
	
	
	esetTmp <- parseProgenesisPeptideCsv(file=progenesisPeptideCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisPeptideCsvFile1))
	sqaGlobal <- safeQuantAnalysis(esetTmp, method="global")
	sqaRt <- safeQuantAnalysis(esetTmp, method="rt")
	stopifnot(sum(sqaRt$pValue < 0.01) > 	sum(sqaGlobal$pValue < 0.01))

	cat("--- testSafequnatAnalysis: PASS ALL TEST --- \n")
	
}


### TEST FUNCTIONS END

### TESTS

testExport() 
testSafequantAnalysis()

#e <- combine(eset[1:10,],eset[20:900,])
#exprs(e)


### TESTS END

