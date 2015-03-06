# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### @TODO 
parseCSV <- function(file=file,expDesign=expDesign){}


################################## LFQ ##################################

###  get column indices of intensity data
#.getProgenesisCsvIntColIndices <- function(file){
#	
#	con <- file(file) 
#	open(con);
#	line1 <- readLines(con, n = 1)
#	close(con)
#	intStartCol <- grep("Normalized abundance",unlist(strsplit(line1,",")))
#	intEndCol <- grep("Raw abundance",unlist(strsplit(line1,",")))-1
#	nbSamples <- intEndCol - intStartCol + 1
#	return((intStartCol:intEndCol)+nbSamples)
#}

##  get column indices of spectral count data
#' @export
.getProgenesisCsvExpressionColIndices <- function(file, method="auc"){
	
	con <- file(file) 
	open(con);
	line1 <- readLines(con, n = 1)
	close(con)
	intStartCol <- grep("Normalized abundance",unlist(strsplit(line1,",")))
	intEndCol <- grep("Raw abundance",unlist(strsplit(line1,",")))-1
	nbSamples <- intEndCol - intStartCol + 1
	
	if(method =="auc"){
		return(return((intStartCol:intEndCol)+nbSamples))
	}else if(method =="spc"){
		
		startCol <- grep("ectral counts",unlist(strsplit(line1,",")))
		endCol <- startCol+nbSamples-1
		return((startCol:endCol))
		
	}else{
		stop("Unknown quantification indices\n")
	}
}


# Get Experimental 	
#' Parse Experimental Design from Progenesis Csv Export
#' @param file path to progenesis csv file
#' @return data.frame describing experimental design
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getExpDesignProgenesisCsv <- function(file){
	
	header <-  read.csv(file,nrows=2, check.names=F)[,.getProgenesisCsvExpressionColIndices(file)]
	header <- data.frame(header) ## in case only one cond one sample
	
	conditionLine <- as.character(unlist(header[1,]))
	currentCond <- as.character(conditionLine[1])
	condition <- c()
	
	for(i in 1:length(conditionLine)){
		
		if(nchar(conditionLine[i]) > 0 ){
			currentCond <- conditionLine[i]
		}
		condition <- c(condition,currentCond) 
		
	}
	
	## set control condition
	isControl <- rep(F,length(condition))
	isControl[condition[1] == condition] <- T
	
	### check if sample names are unique, if not, -> add running number 
	# substitute - -> .
	#sampleNames <- gsub("\\-",".",as.character(as.vector(unlist(header[2,]))))
	sampleNames <- as.character(as.vector(unlist(header[2,])))
	if(T %in% (table(sampleNames) > 1)) sampleNames <- paste(sampleNames,1:length(sampleNames),sep="_") 
	
	expDesign <- data.frame(
			condition=condition
			,isControl=isControl
			,row.names=sampleNames
	)
	return(expDesign)
}

#' Parse Progenesis Protein Csv
#' @param file path to Progenesis Protein csv file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @return ExpressionSet object
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseProgenesisProteinCsv <- function(file=file,expDesign=expDesign, method="auc"){
	
	res <- read.csv(file,skip=2,allowEscapes=T, check.names=F)
	expMatrix <- as.matrix(res[,.getProgenesisCsvExpressionColIndices(file, method=method)])
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))

	row.names(expMatrix) <- res$Accession
	
#	[1] "Accession"                      "Peptide count"                 
#	[3] "Peptides used for quantitation" "Confidence score"              
#	[5] "Anova (p)*"                     "Max fold change"               
#	[7] "Highest mean condition"         "Lowest mean condition"         
#	[9] "Description"                    "A11-03216"                     
	
	featureAnnotations <- data.frame(
			proteinName=res$Accession
			,proteinDescription=res$Description
			,idScore=res[,"Confidence score"]
		#	,nbPeptides=res[,"Peptides used for quantitation"] ### old Progenesis 
			,nbPeptides=res[,which(names(res) == "Unique peptides" | names(res) == "Peptides used for quantitation")[1]] # Progensis QI
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=rep(F,nrow(expMatrix))
			,row.names=res$Accession)
	
	featureAnnotations <- featureAnnotations[!allColNA,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	
	### re-order and exclude channels 
	expMatrix <- as.matrix(expMatrix[!allColNA ,rownames(expDesign)])
	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
}

#' Parse Progenesis Peptide Csv
#' @param file path to Progenesis Peptide csv file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @return ExpressionSet object
#' @export
#' @import affy
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseProgenesisFeatureCsv <- function(file=file,expDesign=getExpDesignProgenesisCsv(file), method="auc"){
	
	### stop if not all samples labelled with a given condition are assigned as control
	# Example:
	#	condition isControl
	#	A11.09066 Condition 1      TRUE
	#	A11.09067 Condition 1     FALSE
	#	A11.09068 Condition 1     FALSE
	#	A11.09070 Condition 2     FALSE
	#	A11.09071 Condition 2     FALSE
	#	A11.09072 Condition 2     FALSE
	if(length(unique(expDesign$condition))  == length(unique(expDesign[!expDesign$isControl,]$condition))){
		stop("Invalid Exp. Design")
	}
	
	# read csv file
	res <- read.csv(file,skip=2,allowEscapes=T,check.names=F)
	expMatrix <- as.matrix(res[,.getProgenesisCsvExpressionColIndices(file, method=method)])
		
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	### Mass filed missing in old Progenesis Feature export
	if(is.null(res$Mass)){res$Mass <- rep(NA,nrow(expMatrix))}
	
	# added
	ptm <- res[,"Variable modifications ([position] description)"]
	nbPtmsPerPeptide <- unlist(lapply(ptm,function(t){
										t <- as.character(t)
										return(sum(unlist(gregexpr("\\|",t)[[1]]) > 0) + (nchar(t) > 0)  )}))
	
	
#	"#"                                              
#	[2] "m/z"                                            
#	[3] "Retention time (min)"                           
#	[4] "Retention time window (min)"                    
#	[5] "Mass"                                           
#	[6] "Charge"                                         
#	[7] "Max fold change"                                
#	[8] "Highest mean condition"                         
#	[9] "Lowest mean condition"                          
#	[10] "Anova"                                          
#	[11] "Maximum CV"                                     
#   [60] "Notes"                                          
#	[61] "Score"                                          
#	[62] "Mass error (u)"                                 
#	[63] "Mass error (ppm)"                               
#	[64] "Protein"                                        
#	[65] "Sequence"                                       
#	[66] "Variable modifications ([position] description)"
#	[67] "Description"  
	
	featureAnnotations <- data.frame(
			proteinName=res$Protein
			,proteinDescription=res$Description
			,peptide=res$Sequence
			,idScore=res$Score
			,mass=res$Mass
			,pMassError=res[,"Mass error (ppm)"]
			,mz=res[,"m/z"]
			,retentionTime=res[,"Retention time (min)"]
			,charge=as.numeric(res$Charge)
			,ptm=ptm
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=rep(F,nrow(expMatrix))
	#		,row.names=res$Accession
			# added
			,nbPtmsPerPeptide = nbPtmsPerPeptide
	)
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	## @TODO what if "X001_Yewtivya" "001_Yewtivya"
	#grepl("X[0-9]",colnames(expMatrix))
	
	### re-order and exclude channels  
	#expMatrix <- expMatrix[,rownames(expDesign)]
	
	#return(createExpressionDataset(expressionMatrix=expMatrix[!allColNA,],expDesign=expDesign,featureAnnotations=featureAnnotations[!allColNA,]))
	
	# discard non peptide annotated rows
	isPep <- nchar(as.character(featureAnnotations$peptide)) > 0 
	featureAnnotations <- data.frame(featureAnnotations)[!allColNA & isPep,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	


	### re-order and exclude channels  
	if(ncol(expMatrix) > 1){
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep ,rownames(expDesign)])
	}else{ # to avoid crash when only one run
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,])
	}

	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
	
}




################################## LFQ END ##############################
	
################################## TMT ##################################

### get line number from which Scaffold file should be read
#' @export
.getSkipLineNb <- function(fileName){
	
	conn<-file(fileName,open="r")
	suppressWarnings(linn<-readLines(conn))
	skip <- 1
	while(regexpr("Bio Sample",as.character(linn[skip])) == -1 ){
		skip <- skip +1 
	}
	close(conn)
	
	return(skip-1)
}

### 6-plex or 10-plex
#' @export
.getNbPlex <- function(fileName){
	
	### parse data
	res <- read.csv(fileName,sep="\t",skip=.getSkipLineNb(fileName),nrows=1)
	return(length(grep("^Normalized",names(res))))
}


### get file type
#' @export
.getFileType <- function(file){
	
	if( sum(grepl("Retention.time..min",names(read.csv(file,skip=2, nrows=0))))  > 0 ){
		return("ProgenesisFeature")
	}else if(sum( grepl("Peptide.count",names(read.csv(file,skip=2, nrows=0)) )) > 0  ){
		return("ProgenesisProtein")
	}else if((grepl("^Identification.Criteria",names(read.csv(file,skip=2, nrows=0))) > 0) &
			(grepl("^Experiment",names(read.csv(file,skip=1, nrows=0))) > 0)){
		return("ScaffoldTMT")
	}else if(F){ 
		return("GenericCSV") 
	}else{
		return("Unknown")	
	}
}


#' Parse scaffold output .xls file (RAW export)
#' @param file path to Scaffold file
#' @param expDesign experimental design data.frame
#' @param keepFirstAcOnly TRUE/FALSE If multiple ACs in Accession.Numbers filed. Then keep the first one only
#' @return ExpressionSet object
#' @export
#' @import affy
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseScaffoldRawFile <- function(fileName, expDesign=expDesign,keepFirstAcOnly=FALSE,isPurityCorrect=T){
	
	### parse data
	res <- read.csv(fileName,sep="\t",skip=.getSkipLineNb(fileName))
	if(keepFirstAcOnly) res$Accession.Numbers <- gsub("\\,.*","",as.character(res$Accession.Numbers))
	res$Accession.Numbers <- gsub("\\[","",as.character(res$Accession.Numbers))
	res$Accession.Numbers <- gsub("\\]","",as.character(res$Accession.Numbers))
	
	###filter out con and decoy
	### HACK, to make sure GFP proteins are not discarded
	res$Accession.Numbers <- gsub(".*P42212.*","sp|P42212|GFP",res$Accession.Numbers)
	

	# CREATE EXPRESSION SET	
	nbPlex <- length(grep("^Normalized",names(res)))
	expMatrix <- as.matrix(res[, 10:(10+(nbPlex-1)) ])
		
	if(isPurityCorrect){
		
#		expMatrix <- log2(expMatrix)
#		expMatrix[!is.finite(expMatrix)] <- NA
		
		expMatrix <- purityCorrectTMT(expMatrix,getImpuritiesMatrix(nbPlex))
		
		
		# expMatrix <- 2^(expMatrix)
	}
	
	### re-order and exclude channels  
#	expMatrix <- expMatrix[,as.numeric(rownames(expDesign))]
#	colnames(expMatrix) <- rownames(expDesign)
		
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	### re-order and exclude channels  
	expMatrix <- as.matrix(expMatrix[!allColNA ,as.numeric(rownames(expDesign))])
	colnames(expMatrix) <- rownames(expDesign)
	
	featureAnnotations <- data.frame(
			proteinName=res$Accession.Numbers
			,proteinDescription=res$Protein.Name 
			,peptide=res$Peptide.Sequence
			#,idScore=res$Score
			#,mass=res$Acquired.Plus.One.Mass-1
			#,pMassError=res$Mass.error..ppm.
			,mz=((as.numeric(as.character(res$Acquired.Plus.One.Mass)) -1)+res$Charge)/res$Charge
			#,retentionTime=res$Retention.time..min.
			,charge=as.numeric(res$Charge)
			,spectrumName = res$Spectrum.Name
			#,ptm=res$Variable.modifications...position..description.
			,isNormAnchor=rep(T,nrow(res))
			,isFiltered=rep(F,nrow(res))
	)
	

	featureAnnotations <- featureAnnotations[!allColNA,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	
	
	
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
}



################################## TMT END ##############################