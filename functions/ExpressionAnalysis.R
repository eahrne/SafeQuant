
### create data frame of per condition  per feature median indices
getMedianIdxPerCond <- function(groupedData){
	
	proteins <- rownames(groupedData[[1]])
	medianIdx <- data.frame(row.names=proteins)
	conditions <- names(groupedData)
	
	for(cond in conditions){
		medianIdx <- data.frame(medianIdx,apply(groupedData[[cond]],1,median))
	}
	names(medianIdx) <- conditions
	rownames(medianIdx) <- proteins
	return(medianIdx)
	
}


### create data frame of per condition per feature summed indices
getSummedIdxPerCond <- function(groupedData){
	
	proteins <- rownames(groupedData[[1]])
	summedIdx <- data.frame(row.names=proteins)
	conditions <- names(groupedData)
	
	for(cond in conditions){
		summedIdx <- data.frame(summedIdx,apply(groupedData[[cond]],1,sum))
	}
	names(summedIdx) <- conditions
	rownames(summedIdx) <- proteins
	return(summedIdx)
	
}

# get identifiers of non-control conditions
getNonControlConditionNbs <- function(groupedData, controlCondition=controlCondition){
	#return(as.numeric(names(groupedData)[!(names(groupedData) == controlCondition)]))
	return(names(groupedData)[!(names(groupedData) == controlCondition)])
}

### get normalization factors with respect to 1st column sample
getNormalizationFactors <- function(raw){
	
	rawDataSums = apply(raw,2, FUN=sum)
	gainFactors = as.numeric(rawDataSums[1]) / as.numeric(rawDataSums)
	
	return(gainFactors)
	
}

### multiply each data column with norm factors
normalizeIntensities <- function(rawData, gainFactors){
	
	normDataColNames = paste("normInt_",names(rawData),sep="")
	normDataMatrix <- matrix(t(apply(rawData, 1, function(.rawData)mapply(gainFactors, .rawData, FUN="*"))),ncol=ncol(rawData),dimnames=list(row.names(rawData),normDataColNames))
	
	return(data.frame( normDataMatrix))
	
}

### create ungroup normalized auc data set 
# 1) get normalizetion factors
# 2) apply normalization factors
# 3) replace missing values with user specified value @TODO derive appropriate value from distribution
createUngroupedNormData <- function(ungroupedRawIntData ,selectedNormAC=selectedNormAC, userSpecifiedMinInt=userSpecifiedMinInt ,verbose=verbose){
	
	
	### select for user specified normalization proteins
	normSelection <- regexpr(selectedNormAC,row.names(ungroupedRawIntData),ignore.case=TRUE) > -1
	
	### if sum == 0, user specified proteins not valid
	if(sum(normSelection) == 0){
		print(paste("Error: Invalid protein selection for normalization",selectedNormAC))
		q(status=-1)
	}	
	
	### get normalization factors
	nomalizationFactors  <- getNormalizationFactors(ungroupedRawIntData[normSelection,])
	
	if(verbose){
		print("Normalization factors")
		print(nomalizationFactors)
	}
	
	### normalize intensities per sample
	ungroupedNormIntData <- normalizeIntensities(ungroupedRawIntData,nomalizationFactors)
	
	if(verbose){
		print(paste("Replace missing values with ",userSpecifiedMinInt))
	}
	
	### replace missing values
	ungroupedNormIntData[ungroupedNormIntData < userSpecifiedMinInt] <- userSpecifiedMinInt
	
	return(ungroupedNormIntData)
	
}


createExpressionDataset <- function(logNormData,expDesign){
	### create expressionDataset
	conditions = c()
	c <-1
	for(cond in names(expDesign)){
		conditions = c(conditions,rep(c,expDesign[,cond]))
		c<-c+1
	}
	
	pData = data.frame(conditions)
	names(pData) = c("condition")
	rownames(pData) = colnames(logNormData)

	metadata = data.frame(labelDescription = "conditon label",row.names= names(pData))
	phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata) 
	annotation = "test"
	
	experimentData = new("MIAME", name = "Erik Ahrne",
			lab = "Proteomics Core Facillity, Biozentrum, Basel", contact = "erik.ahrne@unibas.ch",
			title = "test proteomics bioC", abstract = "An example ExpressionSet",
			url = "www.lab.not.exist", other = list(notes = "Created from csv file"))
	
	eset = new("ExpressionSet", exprs = logNormData,
			phenoData = phenoData, experimentData = experimentData,
			annotation = annotation)
	
	### create expressionDataset end
}

### get differential expression qvalues per condtion, return data frame of qvalues per nonControlConditions
getDEQvaluesPerCondition <- function(groupedNormIntData
		,controlCondition=controlCondition
		,nonControlConditions=nonControlConditions
		,expDesign=expDesign  ){
	expressionMatrix <- log(as.matrix(data.frame(groupedNormIntData[controlCondition],groupedNormIntData[nonControlConditions])))

	expDesignLocal <- data.frame(expDesign[,names(expDesign) %in% controlCondition] ,expDesign[,names(expDesign) %in% nonControlConditions])
	
	names(expDesignLocal) <- c(controlCondition,nonControlConditions)
	
	eset <- createExpressionDataset(expressionMatrix,expDesignLocal)
	e <- exprs(eset)
	f <- factor(as.character(eset$condition))
	design <- model.matrix(~f)
	#### calculate modified t-statistic, empirical Bayes method, Smyth (2004) 
	fit <- eBayes(lmFit(eset,design))
	#### adjust for multiple-testing, Benjamin & Hochberg (1995) method 
	
	if(length(nonControlConditions) >1){
		qvalues <- data.frame(apply(fit$p.value[,2:(length(nonControlConditions)+1)],2,function(pvals){ return(p.adjust(pvals,method="fdr")) }))
		names(qvalues) <- nonControlConditions
		return(qvalues)	
	}else{
		qvalues <- data.frame(p.adjust(fit$p.value[,2],method="fdr"))
		names(qvalues) <- nonControlConditions
		return(qvalues)
	}
	
	
}

### calculate qvalues of differenatial expression for condition pair (control vs. condition)
getPairWiseDEQvaluesPerCondition <- function(groupedNormIntData
		, controlCondition=controlCondition
		, nonControlConditions=nonControlConditions
		, expDesign=expDesign){
	
	proteins <- rownames(groupedNormIntData[[controlCondition]])
	qvalues <- data.frame(row.names=proteins)
	
	### iterate over all nonControlConditions and calc pair wise eBayes
	for(cond in nonControlConditions){
		
		pairedNormIntData <- list()
   		pairedNormIntData[[controlCondition]] <- data.frame(groupedNormIntData[controlCondition])
		pairedNormIntData[[cond]] <- data.frame(groupedNormIntData[cond])
		
		qs <- getDEQvaluesPerCondition(
				pairedNormIntData
				,controlCondition=controlCondition
				,nonControlConditions=cond
				,expDesign=expDesign[,names(expDesign) %in% c(controlCondition,cond)]
    	)
		
		qvalues <- data.frame(qvalues,qs)											
		
	}
	names(qvalues) <- nonControlConditions
	return(qvalues)	
	
}

### create data frame of median ratios per condition   	
getMedianRatiosPerCondition <- function(medianNormIntPerCond,controlCondition=controlCondition,nonControlConditions=nonControlConditions){
	ratios <- data.frame(medianNormIntPerCond[,nonControlConditions] / medianNormIntPerCond[,controlCondition])
	names(ratios) <- nonControlConditions
	return(ratios)
}

#calculate coefficient of variance (sd/mean)
getCV <- function(data){
	return(apply(data,1,sd)/apply(data,1,mean))
}

### calculate coefficient per condition
getCVPerCondition <- function(groupedNormIntData){
	
	conditionNames <- names(groupedNormIntData)
	proteins <- rownames(groupedNormIntData[[conditionNames[1]]])
	cv <- data.frame(row.names = proteins)
	
	for(cond in conditionNames){
		data <- groupedNormIntData[[cond]]
		cv <- data.frame(cv,getCV(data))
	}
	
	names(cv) <- conditionNames
	return(cv)
}

### highlight conditions with missing values
getNAPerCondition <- function(groupedRawIntData, minIntensity=minIntensity){
	
	conditionNames <- names(groupedRawIntData)
	proteins <- rownames(groupedRawIntData[[conditionNames[1]]])
	na <- data.frame(row.names = proteins)
	
	for(cond in conditionNames){
		data <- groupedRawIntData[[cond]]
		condNa <- apply(is.na(data) | (data <= minIntensity),1, sum ) > 0
		na <- data.frame(na,condNa)
	}
	
	names(na) <- conditionNames
	return(na)
}
