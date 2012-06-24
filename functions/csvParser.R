# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


parseProgOutData <- function(progenesisFile, verbose=verbose){
	
	cat("PARSING " ,progenesisFile,"..\n")
	
	if(verbose){
		print(paste("PARSING CSV FILE",progenesisFile))
	}
	
	data <- read.csv(file=progenesisFile,head=TRUE,sep=",", skip=2)
	#data <- read.csv(file=progenesisFile,head=TRUE,sep=",", skip=2, nrows=1328321)

	#@TODO Can speed up the paring by pre-defining column classes	
	#	classTest <- read.csv(file=progenesisFile, header = TRUE, skip=2, nrows = 1000 )
	#	classes <- sapply(classTest, class)
	#	#classes[["X."]] <- "factor"
	#	classes[[1]] <- "factor"
	#	print(names(classes))
	#	print(classes)
	#	print(system.time( data <- read.csv(file=progenesisFile,head=TRUE, skip=2, colClasses="factor" )) )
	#	
	
	if(verbose){
		print("col.names")
		print(names(data))
	}
	
	return(data)
}

### extract all data of a given qunatification index
getUngroupedData <- function(pOutData, columnStart=columnStart, nbSamples=nbSamples){
	
	columnEnd <- columnStart+nbSamples-1
	ungroupedData <- pOutData[,columnStart:columnEnd]
	#rownames(ungroupedData) <- pOutData$Accession
	#colNames <- gsub(".[1-9]$","",paste("raw",names(progOutData)[rawIntColumnStart:rawIntColumnEnd],sep=""))
	colNames <- names(pOutData)[columnStart:columnEnd]
	
	### remove trailing integers
	#colNames <- gsub("\\.[0-9]*$","",perl=T,colNames)
	names(ungroupedData) <- colNames
	
	return(ungroupedData)
	
}


### group data in accordance with experimental desing (sample data per condition)
getCondGroupedSampleData <-function(ungroupedData, expDesign=expDesign){
	
	groupedData <- list()
	condStartCol <- 1
	for(condNb in names(expDesign)){
		condEndCol <- condStartCol+as.numeric(expDesign[,condNb])-1
		#condEndCol <- condStartCol+as.numeric(expDesign[condNb])-1
		groupedData[[condNb]] <- data.frame(ungroupedData[,condStartCol:condEndCol])
		condStartCol <- condEndCol+1
	}
	
	return(groupedData)
}
