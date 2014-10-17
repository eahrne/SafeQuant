
.log2Exprs <- function(eset){
	exprs(eset) <- log2(exprs(eset))
	return(eset)
}

.exp2Exprs <- function(eset){
	exprs(eset) <- 2^(exprs(eset))
	return(eset)
}

.getControlCondition <-function(eset){
	return(unique(as.character(pData(eset)$condition[pData(eset)$isControl]))[1])
}

#' Create ExpressionSet object
#' @param expressionMatrix matrix of expression signals per feature and sample
#' @param expDesign experimental design data.frame
#' @param featureAnnotations data.frame including e.g: Protein Description, Id score etc.
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
createExpressionDataset <- function(expressionMatrix=expressionMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations){
	
	### make sure that only one unique condition is specified as control
	if(length(unique(as.character(expDesign$condition[expDesign$isControl])) ) > 1){
		print(pData(eset))
		stop("Invalid experimental design")	
	}
	
	### phenoData: stores expDesign
	# display with pData(eset): 
	#	condition    type
	#	A_rep_1         A Control
	#	A_rep_2         A Control
	#	B_rep_1         B    Case
	#	B_rep_2         B    Case
	#	C_rep_1         C    Case
	#	C_rep_2         C    Case
	pData <- new("AnnotatedDataFrame", data=expDesign)
	
	
	### featureData:add more data to each feature. E.g: Protein Description, Id score etc.
	
	return(ExpressionSet(assayData = expressionMatrix
					, phenoData =  pData							### yeah this is weid, but gives error if rolling up already 
					# rolled up eset unless colindices are explicitly specified		
					, featureData = new("AnnotatedDataFrame", data= featureAnnotations[,1:ncol(featureAnnotations)])  
			)
	)
}

#' Perform statistical test (mderated t-test), comparing all case to control
#' @param eset ExpressionSet
#' @param adjust TRUE/FALSE adjust for multiple testing using Benjamini & Hochberg  (1995) method 
#' @return ExpressionSet object
#' @export
#' @import limma affy
#' @note  No note
#' @details No details
#' @references Empirical Bayes method, Smyth (2004), \url{http://www.ncbi.nlm.nih.gov/pubmed/16646809} 
#' @seealso \code{\link[limma]{eBayes}}
#' @examples print("No examples")
getAllEBayes <- function(eset=eset, adjust=F){
	
	
	controlCondition <- .getControlCondition(eset)
	caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
	
	pvalues <- data.frame(row.names=featureNames(eset))
	
	for(cC in caseConditions){
		
		esetPair <- eset[,unlist(pData(eset)$condition %in% c(controlCondition,cC))]
		
		f <- factor(as.character(esetPair$condition))
		design <- model.matrix(~f)
		#### calculate modified t-statistic, empirical Bayes method, Smyth (2004) 
		fit <- eBayes(lmFit(esetPair,design))
		
		p <- fit$p.value[,2]
		
		if(adjust){ ### adjust for multiple testing using Benjamini & Hochberg  (1995) method 
			p <- p.adjust(p,method="BH")
		}
		
		pvalues <- cbind(pvalues,p)
	}
	
	names(pvalues) <- caseConditions
	return(pvalues)	
	
}

#' Calculate ratios, comparing all case to control
#' @param eset ExpressionSet 
#' @param method median or mean
#' @return ExpressionSet object
#' @export
#' @import Biobase 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getRatios <- function(eset, method="median", log2=T){
	
	controlCondition <- .getControlCondition(eset)
	caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
	
	ratios <- data.frame(row.names=featureNames(eset))
	
	allConditions <-  pData(eset)$condition
	
	for(cC in caseConditions){
		
		mCase <- exprs(eset)[,allConditions == cC]
		mControl <- exprs(eset)[,allConditions == controlCondition]
		
		if(method == "median"){
			if(sum(allConditions == cC) > 1){
				mCase <- unlist(apply(exprs(eset)[,allConditions == cC],1,median, na.rm=T))
			}
			if(sum(allConditions == controlCondition) > 1){
				mControl <- unlist(apply(exprs(eset)[,allConditions == controlCondition],1,median, na.rm=T))
			}
		}else if(method == "mean"){
			if(sum(allConditions == cC) > 1){
				mCase <- unlist(apply(exprs(eset)[,allConditions == cC],1,mean, na.rm=T))
			}
			if(sum(allConditions == controlCondition) > 1){
				mControl <- unlist(apply(exprs(eset)[,allConditions == controlCondition],1,mean, na.rm=T))
			}
		}else{
			stop("Unknown method ",method)
		}
		ratios <- cbind(ratios,mCase/ mControl)
	}
	
	names(ratios) <- caseConditions
	
	if(log2){
		return(log2(ratios))
	}else{
		return(ratios)
	}
}

#' Calculate Coefficiant of Variance per feature (Relative standard Deviation) per Condition
#' @param eset ExpressionSet
#' @return data.frame of CVs per condition
#' @export
#' @import Biobase 
#' @note  No note
#' @details CV = sd / mean 
#' @references NA
#' @seealso  \code{\link{getCV}}
#' @examples print("No examples")
getAllCV <- function(eset){
	
	allConditions <- unique(pData(eset)$condition)
	cv <- data.frame(row.names=featureNames(eset))
	for(cond in allConditions){
		
		cvTmp <- rep(NA,nrow(cv))
		colMatch <- pData(eset)$condition ==  cond
		
		### if replicates
		if(sum(colMatch) > 1){
			cvTmp <- getCV(exprs(eset)[,colMatch])
		}
		
		cv <- cbind(cv,cvTmp)
	}
	names(cv) <- allConditions
	return(cv)
}

#' Get signal at zscore x (x standard deviations below mean)
#' @param intensities refrence run signals
#' @param percentile baseline value set as specified promille
#' @return baseline value
#' @export
#' @note  No note
#' @references NA
#' @examples print("No examples")
getBaselineIntensity <- function(intensities , promille = 5){
	
	### 
	#intensities <- sample(intensities[!is.na(intensities)],1000,replace=T)
	if((promille < 1) |  (promille > 1001)){
		stop("Invalid percentile: ", promille)
	}
	
	intensities <- intensities[is.finite(intensities)]
	
	suppressWarnings(return(quantile(intensities,probs = seq(0,1,0.001),na.rm=T)[promille+1]))
	
}

#' Calculate Coefficiant of Variance per feature (Relative standard Deviation) 
#' @param data data.frame of replicate signals
#' @return vector of CVs
#' @export
#' @note  No note
#' @details CV = sd / mean 
#' @references NA 
#' @examples print("No examples")
getCV <- function(data){
	return(apply(data,1,sd,na.rm=T)/apply(data,1,mean,na.rm=T))
}

#' Get retentiontime base normalization factors
#' @param eset ExpressionSet
#' @param minFeaturesPerBin  minumum number of features per bin. If nb. features are < minFeaturesPerBin -> include neighbouring bins.
#' @return data.frame normalization factors per retention time bin (minute)
#' @export
#' @import limma affy
#' @note  No note
#' @details No details
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
getRTNormFactors <- function(eset, minFeaturesPerBin=100){
	
	### check if eset contain necessary retentionTime colummn
	if(is.null(fData(eset)$retentionTime)){
		stop("retentionTime not found")
	}
	
	### check if isNormAnchor and isFiltered columns are defiend, if -> get normalization factors from nonFiltered anchor proteins
	if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){
		sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered
		if(sum(sel) == 0){
			return(stop("Invalid Anhcor Protein Selection"))
		}
		eset <- eset[sel,]
	}
	
	### make sure minFeaturesPerBin is not > tot. nb. features
	minFeaturesPerBin <- min(c(minFeaturesPerBin,nrow(eset)))
	
	# get all ratios to sample 1
	# @TODO How to select reference run?
	ratios <- log2(exprs(eset)) - log2(exprs(eset)[,1])
	#ratios <- log2(exprs(eset)) - log2(apply(exprs(eset),1,median,na.rm=T))
	
	### get median ratio per minute bin
	roundedRT <- round(fData(eset)$retentionTime)
	rtBin <- sort(unique(round(fData(eset)$retentionTime)))
	
	# normalization factors per retention time bin bin
	rtNormFactors <- data.frame(row.names=rtBin)
	
	# iterate over samples
	for(i in 1:ncol(ratios)){
		ratio <- ratios[,i]
		rtNFac <- data.frame(row.names=rtBin)
		
		# iterate over all retention time bins
		for(rt in rtBin){
			
			match <- roundedRT %in% rt
			
			### while number of features are < minFeaturesPerBin -> include neighbouring bins
			rtBins <- rt
			while(sum(match) < minFeaturesPerBin ){
				rtBins <- ((min(rtBins)-1):(max(rtBins)+1))
				match <- roundedRT %in% rtBins
			}
			rtNFac[as.character(rt),1] <- median(ratio[match],na.rm=T)
		}
		
		# store rt bin norm factors
		rtNormFactors <- cbind(rtNormFactors,rtNFac)
		
	}
	names(rtNormFactors) <- rownames(pData(eset))
	return(rtNormFactors)
}

#' Normalization data per retention time bin
#' @param eset ExpressionSet
#' @param rtNormFactors  obtained using getRTNormFactors
#' @return data.frame normalization factors per retention time bin (minute)
#' @export ExpressionSet
#' @import limma affy
#' @note  No note
#' @details Normalize for variations in elelctrospray ionization current.
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
rtNormalize <- function(eset,rtNormFactors){
	
	### check if eset contain necessary retentionTime colummn
	if(is.null(fData(eset)$retentionTime)){
		stop("retentionTime not found")
	}
	
	
	roundedRT <-  round(fData(eset)$retentionTime)
	
	for(i in 1:ncol(eset)){
		# normalize 
		exprs(eset)[,i]  <- 2^(log2(exprs(eset)[,i]) - rtNormFactors[as.character(roundedRT),i])
	}
	
	return(eset)
}


#' Get normalization factors. calculated as median signal per run (column) over median of first run. 
#' @param eset ExpressionSet
#' @return vector of normalization factors
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @examples print("No examples")
getGlobalNormFactors <- function(eset, method="median"){
	
	sel <- rep(T,nrow(eset))
	
	### check if isNormAnchor and isFiltered columns are defiend, if -> get normalization factors from nonFiltered anchor proteins
	if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){
		sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered
		if(sum(sel) == 0){
			stop("Error: Invalid Anhcor Protein Selection")
		}
	}
	
	if(method == "median"){
		rawDataIdx = apply(data.frame(exprs(eset)[sel,]),2, FUN=median, na.rm=T)
	}else if(method == "mean"){
		rawDataIdx = apply(data.frame(exprs(eset)[sel,]),2, FUN=mean, na.rm=T)
	}else{
		stop("Error: Unknown Global Normalization Method", method)		
	}
	
	normFactors = as.numeric(rawDataIdx[1]) / as.numeric(rawDataIdx)
	
	return(normFactors)
	
}

#' Normalize, Norm factors calculated as median signal per run (column) over median of first run. 
#' @param eset ExpressionSet
#' @return eset ExpressionSet
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso getGlobalNormFactors
#' @examples print("No examples")
globalNormalize <- function(eset,globalNormFactors){
	
	normData <- data.frame(matrix(t(apply(exprs(eset), 1, function(.rawData)mapply(globalNormFactors, .rawData, FUN="*"))),ncol=ncol(exprs(eset)),dimnames=list(row.names(exprs(eset)),names(exprs(eset)))))
	names(normData) <- sampleNames(eset) 
	
	return(createExpressionDataset(as.matrix(normData),pData(eset),fData(eset)))
	
}


#' Normalize 
#' @param eset ExpressionSet
#' @param	method c("global","rt") 
#' @return eset ExpressionSet
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso getGlobalNormFactors, getRTNormFactors
#' @examples print("No examples")
normalize <- function(eset, method="global"){
	
	esetNorm <- eset
	
	if("global" %in% method){
		
		globalNormFactors <- getGlobalNormFactors(esetNorm)
		
		### add normalization factors to ExpressionSet
		pData(esetNorm) <- cbind(pData(esetNorm),globalNormFactors)
		esetNorm <- globalNormalize(esetNorm,globalNormFactors)
	}
	if("rt" %in% method){
			
			rtNormFactors <- getRTNormFactors(esetNorm, minFeaturesPerBin=100)
			esetNorm <- rtNormalize(eset,rtNormFactors)
			
	}
	

	### @ experimental
#	if("vsn" %in% method){
#		library(vsn)
#		esetNorm <- justvsn(eset)
#		exprs(esetNorm) <- 2^exprs(esetNorm)
#	}
	
	if(identical(esetNorm,eset)){
		warning("No normalization performed")
	}
	
	return(esetNorm)

}


#' Summarize replicate signal per condition (min)
#' @param data data.frame of replicate signals
#' @param method median (default), mean, max, min, sd
#' @return data.frame of per condition signals
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getSignalPerCondition <- function(eset,method="median"){
	
	conditionNames <- levels(pData(eset)$condition)
	perCondSignal <- data.frame(row.names=rownames(eset))
	
	if(method=="median"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,median))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="mean"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,mean))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}
	else if(method=="max"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,max))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="min"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,min))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="sd"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,sd))
			}else{
				perCondSignal <- cbind(perCondSignal,NA)
			}
		}
	}
	else{
		stop("Unknown method: method ")
	}
	names(perCondSignal) <- conditionNames
	
	return(perCondSignal)
}

#' Roll up feature intensites per unique colum combination
#' @param eset ExpressionSet
#' @param featureDataColumnName vector of column names e.g. peptide or proteinName
#' @param method "sum", "mean" or "top3" 
#' @param isProgressBar TRUE/FALSE display progress bar
#' @return ExpressionSet object
#' @details featureDataColumnName = c("peptide","charge","ptm"), method= c("sum"), sums up intensities per peptie modification charge state
#' @export
#' @import Biobase
#' @note  No note
#' @references No references
#' @seealso \code{\link{topX}}
#' @examples print("No examples")
rollUp <- function(eset=eset,featureDataColumnName= c("peptide"), method=c("sum"),isProgressBar=F ){
	
	# find columns matching featureDataColumnName
	selectedColumns <- names(fData(eset)) %in% featureDataColumnName 
	if(sum(selectedColumns) == 0){
		stop("Unknown featureDataColumnName ",featureDataColumnName,"\n")
	}
	
	### create index tags by concatinating selected coulm entries
	allIndexTags <- as.vector(unlist(apply(data.frame(fData(eset)[,selectedColumns]),1,function(t){
								return(paste(as.vector(unlist(t)),collapse="_"))
							})))
	
	uniqueIndexTags <- unique(allIndexTags)
	nbUniqueIndexTags <- length(uniqueIndexTags)
	
	### data frame of rolled up feature data and assay data
	rolledFData <- data.frame( matrix(ncol=ncol(fData(eset)),nrow=0))
	names(rolledFData) <- names(fData(eset))
	rolledAssayData <- data.frame( matrix(ncol=ncol(eset),nrow=nbUniqueIndexTags) ,row.names=uniqueIndexTags )
	names(rolledAssayData) <- colnames(eset)
	
	### additional columns to be added to rolledFData
	nbRolledFeatures <- vector(length=nbUniqueIndexTags)
	allChargeStates <- rep(NA,nbUniqueIndexTags)
	allPtms <- rep(NA,nbUniqueIndexTags)
	allSpectrumNames <- rep(NA,nbUniqueIndexTags)
	
	### progress bar
	pbSum <- txtProgressBar(min = 0, max = nbUniqueIndexTags, style = 3)
	
	for(i in 1:nbUniqueIndexTags){
		
		### increment progress bar
		if(isProgressBar) setTxtProgressBar(pbSum, i)
		
		matchIndices <- which(allIndexTags %in% uniqueIndexTags[i])
		
		### assay data
		aDataSubset <- exprs(eset)[matchIndices,]
		
		# feature data
		fDataSubset <- fData(eset)[matchIndices,]
		nbFeat <- length(matchIndices)
		
		rolledAssayDataSubset <- aDataSubset
		rolledFDataSubset <- fDataSubset
		
		### more than one matching rows
		if(nbFeat>1){
			
			### roll up method
			### sum
			if(method =="sum"){
				rolledAssayDataSubset <- apply(aDataSubset,2,sum,na.rm=T)
			}else if(method=="mean") {
				### mean
				rolledAssayDataSubset <- apply(aDataSubset,2,mean,na.rm=T)
			}else if(method== "top3"){
				### top3
				rolledAssayDataSubset <- getTopX(aDataSubset,topX=3)
			}else if(method== "top1"){
				### top1 TEST OPTION
				rolledAssayDataSubset <- getTopX(aDataSubset,topX=1)
			}else{
				stop("Unknown roll up method ", method, "\n")
			}
			
			#@TODO How to best rollUp id scores? Currently taking highest score per roll-up group
			# keep feature data of best feature (order by idScore)
			# if no idScore keep first one
			bestIndex <- 1
			if(!is.null(fDataSubset$idScore)){
				bestIndex <- order(fDataSubset$idScore, decreasing = T)[1]
			}
			rolledFDataSubset <- fDataSubset[bestIndex,]
			
		}	
		
		# add rolled data
		rolledAssayData[i,] <- rolledAssayDataSubset
		rolledFData <- rbind(rolledFData,rolledFDataSubset)
		nbRolledFeatures[i] <- nbFeat
		if(!is.null(fDataSubset$charge)){ ### all charge state tag e.g. "2:3:4"
			allChargeStates[i] <- paste(sort(unique(fDataSubset$charge)),collapse=":")
		}
		if(!is.null(fDataSubset$ptm)){ ### all ptms tag e.g. "[5] Phospho (ST):[3] Phospho (ST)|[5] Phospho (ST)"  
			allPtms[i] <- paste(unique(fDataSubset$ptm),collapse=":")
		}
		if(!is.null(fDataSubset$spectrumName)){ ### all spectrum name tag e.g. "07_07Da.55018.55018.3:007_07Da.55013.55013.2:008_07Da.43303.43303.2:006_07Da.54094.5"  
			allSpectrumNames[i] <- paste(unique(fDataSubset$spectrumName),collapse=":")
		}
		
	}
	
	# close progress bar
	setTxtProgressBar(pbSum, i)
	close(pbSum)
	
	### format rolledFData and add rolled up summary columns
	rownames(rolledFData) <- rownames(rolledAssayData)
	rolledAssayData <- as.matrix(rolledAssayData)
	rolledAssayData[rolledAssayData == 0 ] <- NA
	#names(rolledFData) <- paste("best",names(fData(eset)),sep="_")
	names(rolledFData) <- names(fData(eset))
	rolledFData <- cbind(rolledFData,nbRolledFeatures,allChargeStates,allPtms,allSpectrumNames)
	# reset anchor
	#if(!is.null(rolledFData$best_isNormAnchor)){
	if(!is.null(rolledFData$isNormAnchor)){
		#	names(rolledFData)[names(rolledFData) == "best_isNormAnchor"] <- "isNormAnchor"
		fDataSubset$isNormAnchor <- T
	}
	# reset filter
	#if(!is.null(rolledFData$best_isFiltered)){
	if(!is.null(rolledFData$isFiltered)){
		#names(rolledFData)[names(rolledFData) == "best_isFiltered"] <- "isFiltered"
		fDataSubset$isFiltered <- F
	}
	
	return(createExpressionDataset(expressionMatrix=rolledAssayData,expDesign=pData(eset),featureAnnotations=rolledFData))
}

#' Calculate Mean of X most intense features
#' @param entryData data.frame listing feature intensities of one entry 
#' Typically rows corresponds to Peptide entries of one protein
#' @param topX best X flyers
#' @return vector of topX intensities per column (sample) 
#' @export
#' @import  
#' @note  No note
#' @details No details
#' @references Absolute quantification of proteins by LCMSE: A virtue of parallel MS acquisition, Silva (2006), \url{http://www.ncbi.nlm.nih.gov/pubmed/16219938}, 
#' Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Ahrne (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23794183}
#' @examples print("No examples")
getTopX <- function(entryData,topX=3){
	
	# if number of rows are fewer than topX  
	topX <- min(c(topX,nrow(entryData)))
	
	### get row order based on decreasing row sum
	if( !is.null(ncol(entryData))  && (ncol(entryData) > 1)  ){
		o <- order(apply(entryData,1,sum,na.rm=T),decreasing=T)[1:topX]
		
		if(topX==1){
			return(entryData[o,])
		}
		return(apply(entryData[o,],2,mean,na.rm=T))
	}
	
	# only one column
	o <- order(entryData,decreasing=T)[1:topX]
	return(mean(entryData[o],na.rm=T))	
	
}

#' Calculate intensity-based absolute-protein-quantification (iBAQ) metric per protein 
#' @param eset protein level ExpressionSet
#' @param list protein sequneces
#' @param peptideLength peptide length interval (to get number of peptides used for normalization)
#' @param nbMiscleavages number of mis-cleavages allowed when digesting protein sequneces in silico (to get number of peptides used for normalization)
#' @param proteaseRegExp protease Reg Exp cleavage rule
#' @return ExpressionSet
#' @export
#' @import  
#' @note  No note
#' @details No details
#' @references Global quantification of mammalian gene expression control, Schwanhausser (2011), \url{http://www.ncbi.nlm.nih.gov/pubmed/21593866}, 
#' Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Ahrne (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23794183}
#' @examples print("No examples")
getIBAQEset <- function(eset
		, proteinDB=NA
		, peptideLength = c(5,36)
		, nbMiscleavages = 0
		, proteaseRegExp=.getProteaseRegExp("trypsin")
){
		
	# get number of detectable peptides per protein
	nbPeptides <- vector(length=nrow(eset))
	i <- 0
	for(i in 1:nrow(eset)){
		
		ac <- as.character(fData(eset)$proteinName[i])
		
		nbPep <- NA
		if( !is.null(proteinDB[[ac]]) ){
			nbPep <- getNbDetectablePeptides(getPeptides(proteinDB[[ac]],proteaseRegExp=proteaseRegExp,nbMiscleavages=nbMiscleavages),peptideLength=peptideLength)
		}else{
			warning(ac," NOT FOUND IN PROTEIN DATABASE")
		}
		nbPeptides[i] <- nbPep
	}
	
	### create new "absquant" eset
	esetIBAQ <- eset
	# scale protein intensity by number of detectable peptides 
	exprs(esetIBAQ) <- exprs(eset) / nbPeptides
	
	### store normalization factors
	fData(esetIBAQ)$nbTrypticPeptides <- nbPeptides
	
	return(esetIBAQ)
	
}

### require colnames "signal", "cpc"
#' Leave-One-Out Cross Validate Qunatification Model
#' @param data.frame  of two columns 1) "signal" - ms metric 2) "cpc" absolute quantity
#' @return data.frame of fold errors per (left-out) protein 
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso NA
#' @examples print("No examples")
getLoocvFoldError <- function(df){
	
	foldError <- vector(length=nrow(df))
	
	for (i in 1:nrow(df)) {
		fit <- lm(cpc ~ signal, data=df[-i,] )
		foldError[i] <- 10^predict(fit,newdata=df[i,]) / 10^df[i,]$cpc
		
	}
	
	### if foldErro < 1, take neg. of inverse	
	foldError[foldError < 1] <- -1/foldError[foldError < 1]
	
	### return fold error per protein
	foldError <- data.frame(foldError,row.names=rownames(df))
	
	return(foldError)
	
}


#' Set value to NA if it deviatves with more than 1.5 * IQR from lower/upper quantile
#' @param vector numeric
#' @return vector numeric
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso NA
#' @examples print("No examples")
removeOutliers <- function(x, na.rm = TRUE, ...){
	qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
	H <- 1.5 * IQR(x, na.rm = na.rm)
	y <- x
	y[x < (qnt[1] - H)] <- NA
	y[x > (qnt[2] + H)] <- NA
	
	return(y)
}
