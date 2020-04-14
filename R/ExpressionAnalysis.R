
## get index of max (created for data.table)
#' get index of max in vecotor of numeric values
#' @param v vector
#' @export
getMaxIndex <- function(v){
  if(all(is.na(v))) return(1)
  return(which(max(v,na.rm=T) == v)[1])
}

## return first entry per column
#' @export
.getFirstEntry <- function(x){
  return(x[1])
}

#' @export
.log2Exprs <- function(eset){
  exprs(eset) <- log2(exprs(eset))
  return(eset)
}

#' @export
.exp2Exprs <- function(eset){
  exprs(eset) <- 2^(exprs(eset))
  return(eset)
}

#' @export
.getControlCondition <-function(eset){
  return(unique(as.character(pData(eset)$condition[pData(eset)$isControl]))[1])
}

#' Create ExpressionSet object
#' @param expressionMatrix matrix of expression signals per feature and sample
#' @param expDesign experimental design data.frame
#' @param featureAnnotations data.frame including e.g: Protein Description, Id score etc.
#' @return ExpressionSet object
#' @import Biobase
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
createExpressionDataset <- function(expressionMatrix=expressionMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations){

  ### make sure that only one unique condition is specified as control
  if(length(unique(as.character(expDesign$condition[expDesign$isControl])) ) > 1){
    stop("ERROR: createExpressionDataset, Invalid experimental design")
  }

  ### phenoData: stores expDesign
  # display with pData(eset):
  #	condition    isControl
  #	A_rep_1         A T
  #	A_rep_2         A T
  #	B_rep_1         B    F
  #	B_rep_2         B    F
  #	C_rep_1         C    F
  #	C_rep_2         C    F
  #pData <- new("AnnotatedDataFrame", data=expDesign)
  expDesign$condition <- as.factor(gsub(" ","",	expDesign$condition ))
  pData <- AnnotatedDataFrame(data=expDesign)

  ### featureData:add more data to each feature. E.g: Protein Description, Id score etc.

  return(ExpressionSet(assayData = expressionMatrix
                       , phenoData =  pData							### yeah this is weird, but gives error if rolling up already
                       # rolled up eset unless colindices are explicitly specified
                       #, featureData = new("AnnotatedDataFrame", data= featureAnnotations[,1:ncol(featureAnnotations)])
                       , featureData = AnnotatedDataFrame(data= featureAnnotations[,1:ncol(featureAnnotations)])
  )
  )
}


#' Perform statistical test (mderated t-test), comparing all case to control
#' @param eset ExpressionSet
#' @param adjust TRUE/FALSE adjust for multiple testing using Benjamini & Hochberg  (1995) method
#' @param log T/F log-transform expression values
#' @param method c("all","pairwise")
#' @param adjustFilter matrix T/F do not adjust for multiple testing
#' @return data.frame of pvalues per condition comparison
#' @export
#' @import limma Biobase
#' @importFrom stats model.matrix p.adjust
#' @note  No note
#' @details No details
#' @references Empirical Bayes method, Smyth (2004), \url{http://www.ncbi.nlm.nih.gov/pubmed/16646809}
#' @seealso \code{\link[limma]{eBayes}}
#' @examples print("No examples")
getAllEBayes <- function(eset=eset, 
                         adjust=F, 
                         log=T, 
                         method="pairwise",
                         adjustFilter=matrix(F,nrow=nrow(eset),ncol=length(levels(pData(eset)$condition))-1  )
                         ,...){

  #### calculate moderated t-statistic, empirical Bayes method, Smyth (2004)

  # Question: limma multiple groups comparison produces different pvalue comparing with two group comparsion
  # https://support.bioconductor.org/p/44216/
  # https://support.bioconductor.org/p/60556/
  uniqueConditions <- unique(pData(eset)$condition)
  nbUniqueConditions <- length(uniqueConditions)
  controlCondition <- .getControlCondition(eset)
  caseConditions <- setdiff(uniqueConditions ,controlCondition)

  # if no condition has replicates return NA's
  if(max(table(pData(eset)$condition)) == 1){
    return(data.frame( matrix(nrow=nrow(eset), ncol = length(caseConditions), dimnames = list(rownames(eset),caseConditions))))
  }

  if(log)(exprs(eset) <- log2(exprs(eset)))

  #	#### NON-PAIR-WISE COMPARISONS
  if("all" %in% method){

    # add subject term to allow for paired t-statistic
    if("subject" %in% names(pData(eset))){
      design <- model.matrix(~0+condition + subject, data=pData(eset))
    }else{
      # no intercept
      design <- model.matrix(~0+condition, data=pData(eset))
    }
    colnames(design) <- gsub("^condition","",colnames(design))

    fit <- lmFit(eset,design)
    ###  create contract matrix describing desired condition comparisons
    # e.g. B-A, C-A
    #		contrastMatrix
    #		#    B  C
    #		# A -1 -1
    #		# B  1  0
    #		# C  0  1
    contrastMatrix <- makeContrasts(contrasts=paste(caseConditions,"-", controlCondition),levels=design)
    colnames(contrastMatrix) <- caseConditions

    # fit contrasts coefficients
    fitContrasts <- eBayes(contrasts.fit(fit,contrastMatrix),
                           ...
                           )
    pvalues <- data.frame(fitContrasts$p.value[,caseConditions])
    names(pvalues) <- caseConditions

    #print(head(round(fitContrasts$coefficients[,caseConditions],3) == round(compRatios[,caseConditions],3)))
    #print(cbind(fitContrasts$coefficients[,caseConditions[1]],compRatios[,caseConditions[1]]))

  }else{ 	#### PAIR-WISE COMPARISONS "pairwise" %in% method

    ## REASONING
    # General case (linear model, anova ..)
    # Adding additional groups alters all Std. Errors. However if equal variance (homoscedastic) the std.error and p-values should not increase.
    #   -> We are increasing the number of samples used to estimate the standard error of the condition coefficients
    #	-> If equal variance the power increases
    #
    # limma
    # The eBayes() is also affected by adding more groups as more data is used to estimated the variance prior (derived from data across genes).

    pvalues <- data.frame(row.names=featureNames(eset))

    for(cC in caseConditions){

      ### at least one replicate of one condition requires
      selCol <- unlist(pData(eset)$condition %in% c(controlCondition,cC))
      if(sum(selCol) > 2 ){
        esetPair <- eset[,selCol]
        # add subject term to allow for paired t-statistic
        if("subject" %in% names(pData(eset))){
          fit <-eBayes(lmFit(esetPair, model.matrix(~factor(esetPair$condition) + subject, data=pData(esetPair))),
                       ...)
        }else{
          fit <-eBayes(lmFit(esetPair, model.matrix(~factor(esetPair$condition), data=pData(esetPair)),
                       ),...)
        }
        p <- fit$p.value[,2]

      }else{
        p <- rep(NA,nrow(eset))
      }
      pvalues <- cbind(pvalues,p)
    }
    names(pvalues) <- caseConditions
  }

  if(adjust){ ### adjust for multiple testing using Benjamini & Hochberg (1995) method
    for(i in 1:ncol(pvalues)){

      pvaluesTmp <-  pvalues[,i]
      pvaluesTmp[adjustFilter[,i]] <- NA
      pvalues[,i] <-p.adjust(pvaluesTmp,method="BH")
    }
  }
  return(pvalues)
}

#' Calculate ratios, comparing all case to control
#' @param eset ExpressionSet
#' @param method median, mean, paired
#' @param log2 transform
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getRatios <- function(eset, method="median", log2=T){

  if(sum(c("median","mean","paired") %in% method)  == 0 ){
    stop("Unknown method ",method)
  }
  if(("paired" %in% method) & !"subject" %in% names(pData(eset))){
    stop("Invalid method ",method)
  }

  controlCondition <- .getControlCondition(eset)
  caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
  ratios <- data.frame(row.names=featureNames(eset))
  eCtrl <- subset(exprs(eset),select=which(pData(eset)$condition == controlCondition))

  for(cC in caseConditions){
    eCase <- subset(exprs(eset),select=which(pData(eset)$condition == cC))

    if("subject" %in% names(pData(eset))){

      if(log2){
        rTmp <- log2(eCase) - log2(eCtrl)
      }else{
        rTmp <- eCase /	eCtrl
      }

      if("paired" %in% method){ ### return paired all paired ratios instead of mean/median per condition
        ratios <- cbind(ratios,rTmp)
      }else{
        ratios <- cbind(ratios,apply(rTmp,1,median))
      }
    }else{ # non-paired design
      if(log2){
        ratios <- cbind(ratios,apply(log2(eCase),1,method, na.rm=T) - apply(log2(eCtrl),1,method, na.rm=T))
      }else{
        ratios <- cbind(ratios,apply(eCase,1,method, na.rm=T) /	apply(eCtrl,1,method, na.rm=T))
      }
    }
  }

  if("paired" %in% method){
    names(ratios) <- rownames(pData(eset))[!pData(eset)$isControl]
  }else{
    names(ratios) <- caseConditions
  }

  return(ratios)
}

#' Calculate Coefficiant of Variance per feature (Relative standard Deviation) per Condition
#' @param eset ExpressionSet
#' @return data.frame of CVs per condition
#' @import Biobase
#' @export
#' @note  No note
#' @details CV = sd / mean
#' @references NA
#' @seealso  \code{\link{getCV}}
#' @examples print("No examples")
getAllCV <- function(eset){


  signal <- exprs(eset)
  conditions <- unique(pData(eset)$condition)
  expDesign <- pData(eset)
  cv <- data.frame(row.names=featureNames(eset))

  #if paired design
  if("subject" %in% names(pData(eset))){
    #return(data.frame((matrix(nrow = nrow(eset), ncol = length(allConditions), dimnames=list(rownames(eset), allConditions) ))))

    controlCondition <- .getControlCondition(eset)
    caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
    conditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)

    # discard control from expDesign
    #expDesign <- subset(expDesign, condition %in% conditions) # not accepted by CRAN
    expDesign <- expDesign[expDesign$condition %in% conditions,]

    # add all NAs for control condition
    cv <- cbind(cv,rep(NA,nrow(cv)))
    names(cv) <- controlCondition

    signal <- getRatios(eset,method="paired",log2=F)

  }

  for(cond in conditions){

    cvTmp <- rep(NA,nrow(cv))
    colMatch <- expDesign$condition ==  cond

    ### if replicates
    if(sum(colMatch) > 1){
      cvTmp <- getCV(subset(signal, select=colMatch))
    }

    cv <- cbind(cv,cvTmp)
    names(cv)[ncol(cv)] <- cond
  }

  return(cv)
}

#' Get signal at zscore x (x standard deviations below mean)
#' @param intensities refrence run signals
#' @param promille baseline value set as specified promille
#' @return baseline value
#' @export
#' @importFrom stats quantile
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
#' @importFrom stats sd
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
#' @import limma Biobase
#' @export
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
  #	if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){
  #		sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered
  #		if(sum(sel) == 0){
  #			return(stop("Invalid Anchor Protein Selection"))
  #		}
  #		eset <- eset[sel,]
  #	}
  # No big difference if using all features. Avoids down stream bug when applying norm factors

  ### make sure minFeaturesPerBin is not > tot. nb. features
  minFeaturesPerBin <- min(c(minFeaturesPerBin,nrow(eset)))

  # get all ratios to sample 1
  # @TODO How to select reference run?
  #ratios <- log2(exprs(eset)) - log2(exprs(eset)[,1])
  ratios <- log2(exprs(eset)) - apply(log2(exprs(eset)),1,mean,na.rm=T)

  ### get median ratio per minute bin
  roundedRT <- round(fData(eset)$retentionTime)
  rtBin <- sort(unique(round(fData(eset)$retentionTime)))
  rtBin <-  rtBin[!is.na(rtBin)]

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
      # median or sum?
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
#' @export
#' @import limma Biobase
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

  # FIX
  if(sum(!(as.character(roundedRT)[!is.na(roundedRT)] %in% row.names(rtNormFactors))) > 0){
    stop("ERROR: Missing R.T. Norm Factors.")
  }

  for(i in 1:ncol(eset)){
    # normalize

    nF <- rtNormFactors[as.character(roundedRT),i]
    nF[is.na(nF)] <- 0

    exprs(eset)[,i]  <- 2^(log2(exprs(eset)[,i]) - nF)
  }

  return(eset)
}

#' Get normalization factors. calculated as summed/median signal per run (column) over summed/median of first run.
#' @param eset ExpressionSet
#' @param method c("sum","median", "keepCDiff")
#' @param keepCondDiff logical default False, only normalize within conditions  
#' @return vector of normalization factors
#' @export
#' @note  No note
#' @details method - 'keepCDiff' logical default False, only normalize within conditions  
#' @keywords normalization
#' @references NA
#' @examples print("No examples")
getGlobalNormFactors <- function(eset, method="median" ){

  keepCondDiff = ifelse("keepCDiff" %in% method, T,F)  
  
  if("sum" %in% method){
    method <- "sum"
  }else{
    method <- "median"
  }

  sel <- rep(T,nrow(eset))

  ### check if isNormAnchor and isFiltered columns are defined, if -> get normalization factors from nonFiltered anchor proteins
  if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){

    ### only use feature qunatified in all runs for normalization
    isAllSignal <- as.vector(apply(is.finite(exprs(eset)),1,sum, na.rm=T) == ncol(eset))

    #isCalMix <- fData(eset)$proteinName %in% names(CALIBMIXRATIOS)
    #sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered & isAllSignal & !isCalMix

    sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered & isAllSignal

    if(sum(sel) == 0){
      #stop("Error: getGlobalNormFactors -> Invalid Anchor Protein Selection ")
      cat("WARNING: getGlobalNormFactors -> No proteins matching Anchor Protein Selection \n ")
      return(rep(1,ncol(eset)))
    }
  }

  ### equalize all runs 
  #eD = data.frame(exprs(eset)[sel,])
  eD = exprs(eset)[sel,]
  # if only one row (1 selected protein | peptide for norm)
  if(class(eD) == 'numeric'){
    rawDataIdx = eD
  }else{
    rawDataIdx = apply(eD,2, FUN=method, na.rm=T)
  }
  normFactors = as.numeric(rawDataIdx[1]) / as.numeric(rawDataIdx)
  
  # keep differences between conditions
  if(keepCondDiff){
    # create rawDataIdx Tibble
    if(keepCondDiff){
      rawDataIdxTbl = tibble(cond=eset$condition,rawDataIdx)
      # equalize per condition
      if("median" %in% method){
        condNormFactors = group_by(rawDataIdxTbl, cond) %>% summarise(condNF=median(rawDataIdx)) 
      }else{
        condNormFactors = group_by(rawDataIdxTbl, cond) %>% summarise(condNF=sum(rawDataIdx)) 
      }
      condNormFactors$condNF = condNormFactors$condNF / condNormFactors$condNF[1]
      
      # adjust normFactors by condition norm factors to preserve orignial diff btwn condition median/sums
      for(cond in unique(eset$condition)){
        selCond = eset$condition %in% cond
        normFactors[selCond ] = normFactors[selCond ] * condNormFactors[  condNormFactors$cond == cond, ]$condNF
      }
    }
  }

  return(normFactors)

}

#' Normalize, Norm factors calculated as median signal per run (column) over median of first run.
#' @param eset ExpressionSet
#' @param globalNormFactors globalNormFactors
#' @return eset ExpressionSet
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA
#' @seealso getGlobalNormFactors
#' @examples print("No examples")
globalNormalize = function(eset,globalNormFactors){
  normData =  sweep( exprs(eset), 2, globalNormFactors, FUN="*")
  colnames(normData) <- sampleNames(eset)
  return(createExpressionDataset(normData,pData(eset),fData(eset)))
}

# globalNormalize <- function(eset,globalNormFactors){
#   normData <- data.frame(matrix(t(apply(exprs(eset), 1, function(.rawData)mapply(globalNormFactors, .rawData, FUN="*"))),ncol=ncol(exprs(eset)),dimnames=list(row.names(exprs(eset)),names(exprs(eset)))))
#   names(normData) <- sampleNames(eset)
#   return(createExpressionDataset(as.matrix(normData),pData(eset),fData(eset)))
#
# }

#' Normalize
#' @param eset ExpressionSet
#' @param	method c("global","rt","quantile","keepCDiff)
#' @return eset ExpressionSet
#' @export
#' @import limma
#' @note  No note
#' @details method - 'keepCDiff' logical default False, only normalize within conditions.  Works in cobinaiton witb 'global'
#' @keywords normalization
#' @references NA
#' @seealso getGlobalNormFactors, getRTNormFactors
#' @examples print("No examples")
sqNormalize <- function(eset, method="global"){

  esetNorm <- eset

  if("global" %in% method){

    globalNormFactors <- getGlobalNormFactors(esetNorm,method=method)
    ### add normalization factors to ExpressionSet
    pData(esetNorm)$globalNormFactors <- globalNormFactors
    esetNorm <- globalNormalize(esetNorm,globalNormFactors)
  }
  if("rt" %in% method){

    rtNormFactors <- getRTNormFactors(esetNorm, minFeaturesPerBin=100)
    esetNorm <- rtNormalize(esetNorm,rtNormFactors)

  }
  if("quantile" %in% method){
    # limma
    exprs(esetNorm) <- normalizeQuantiles(exprs(esetNorm))
  }

  ### @ experimental
  if("vsn" %in% method){
  		library(vsn)
  		esetNorm <- justvsn(eset)
  		exprs(esetNorm) <- 2^exprs(esetNorm)
  }

  #	if(identical(esetNorm,eset)){
  #		warning("No normalization performed")
  #		print( apply(exprs(esetNorm),2,median,na.rm=T))
  #	}

  return(esetNorm)

}

#' Summarize replicate signal per condition (min)
#' @param eset ExpressionSet
#' @param method median (default), mean, max, min, sd, sum
#' @return data.frame of per condition signals
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getSignalPerCondition <- function(eset,method="median"){

  conditionNames <- unique(as.character(pData(eset)$condition))
  perCondSignal <- data.frame(row.names=rownames(eset))

  if(sum(c("median","mean","mean","max", "min","sd","sum") %in% method)  == 0 ){
    stop("Unknown method ",method)
  }

  for(cond in conditionNames){
    condMatch <-  cond== pData(eset)$condition
    perCondSignal <- cbind(perCondSignal,apply(subset(exprs(eset),select=which(condMatch)),1,method,na.rm=T))
  }
  names(perCondSignal) <- conditionNames

  return(perCondSignal)
}


#	}
#
#
#
#	else if(method=="medianOld"){
#		for(cond in conditionNames){
#			condMatch <-  cond== pData(eset)$condition
#			if(sum(condMatch) > 1){ # @TODO replace by subset function
#				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,median, na.rm=T))
#			}else{
#				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
#			}
#		}
#	}
#
#
#
#
#	else if(method=="mean"){
#		for(cond in conditionNames){
#			condMatch <-  cond== pData(eset)$condition
#			if(sum(condMatch) > 1){
#				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,mean, na.rm=T))
#			}else{
#				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
#			}
#		}
#	}
#	else if(method=="max"){
#		for(cond in conditionNames){
#			condMatch <-  cond== pData(eset)$condition
#			if(sum(condMatch) > 1){
#				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,max, na.rm=T))
#			}else{
#				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
#			}
#		}
#	}else if(method=="min"){
#		for(cond in conditionNames){
#			condMatch <-  cond== pData(eset)$condition
#			if(sum(condMatch) > 1){
#				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,min, na.rm=T))
#			}else{
#				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
#			}
#		}
#	}else if(method=="sd"){
#		for(cond in conditionNames){
#			condMatch <-  cond== pData(eset)$condition
#			if(sum(condMatch) > 1){
#				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,sd, na.rm=T))
#			}else{
#				perCondSignal <- cbind(perCondSignal,NA)
#			}
#		}
#	}
#	else{
#		stop("Unknown method: method ")
#	}




#' Calculate Mean of X most intense features
#' @param entryData data.frame listing feature intensities of one entry. Typically rows corresponds to Peptide entries of one protein
#' @param topX best X flyers
#' @return vector of topX intensities per column (sample)
#' @export
#' @note  No note
#' @details No details
#' @references Absolute quantification of proteins by LCMSE: A virtue of parallel MS acquisition, Silva (2006), \url{http://www.ncbi.nlm.nih.gov/pubmed/16219938},  Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Ahrne (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23794183}
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
#' @param proteinDB list protein sequneces
#' @param peptideLength peptide length interval (to get number of peptides used for normalization)
#' @param nbMiscleavages number of mis-cleavages allowed when digesting protein sequneces in silico (to get number of peptides used for normalization)
#' @param proteaseRegExp protease Reg Exp cleavage rule
#' @return ExpressionSet
#' @export
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

    ### keep first listed entry sp|P04049|RAF1_HUMAN;sp|P15056|BRAF_HUMAN  -> sp|P04049|RAF1_HUMAN
    #ac <- gsub("\\;.*","",ac)

    nbPep <- NA
    if( !is.null(proteinDB[[ac]]) ){
      nbPep <- getNbDetectablePeptides(getPeptides(proteinDB[[ac]],proteaseRegExp=proteaseRegExp,nbMiscleavages=nbMiscleavages),peptideLength=peptideLength)
    }else{
      warning("WARN: ",ac," NOT FOUND IN PROTEIN DATABASE")
      #cat("WARN: ",ac," NOT FOUND IN PROTEIN DATABASE\n")
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
#' @param df data.frame  of two columns 1) "signal" - ms metric 2) "cpc" absolute quantity
#' @return data.frame of fold errors per (left-out) protein
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA
#' @seealso NA
#' @examples print("No examples")
getLoocvFoldError <- function(df){

  ok <- is.finite(df[,1]) & is.finite(df[,2])
  df <- df[ok,]

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
#' @param x vector numeric
#' @param na.rm logical indicating whether missing values should be removed.
#' @param ... qunatile args
#' @export
#' @importFrom stats quantile IQR
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


#' Roll up feature intensites per unique colum combination
#' @param eset ExpressionSet
#' @param featureDataColumnName vector of column names e.g. peptide or proteinName
#' @param method "sum", "mean" or "top3"
#' @return ExpressionSet object
#' @export
#' @details featureDataColumnName = c("peptide","charge","ptm"), method= c("sum"), sums up intensities per peptie modification charge state
#' @import Biobase data.table
#' @note  No note
#' @references No references
#' @examples print("No examples")
rollUpDT <- function(eset, method = "sum", 	featureDataColumnName =  c("proteinName") ){

  # HACK to please CRAN CHECK "rollUp: no visible binding for global variable "idScore""
  idx <- idScore <- V1 <- allAccessions <- NULL

  ### apply filter
  eset <- eset[!fData(eset)$isFiltered,]

  selectedColumns <- names(fData(eset)) %in% featureDataColumnName
  allIndexTags <- as.vector(unlist(apply(data.frame(fData(eset)[,selectedColumns]),1,function(t){
    return(paste(as.vector(unlist(t)),collapse="_"))
  })))

  DT <- data.table(
    idx = allIndexTags
    ,exprs(eset)
    ,fData(eset)
  )
  setkey(DT,idx)

  if(method =="sum"){
    rDT <- DT[, lapply(.SD, sum, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
    rolledAssayData <- data.frame(rDT,row.names=1, check.names=F)
  }else if(method =="median"){
    rDT <- DT[, lapply(.SD, median, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
    rolledAssayData <- data.frame(rDT,row.names=1, check.names=F)

  }else if(method =="mean"){
    rDT <- DT[, lapply(.SD, mean, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
    rolledAssayData <- data.frame(rDT,row.names=1, check.names=F)
  }else if(method =="top3"){
    rDT <- DT[, lapply(.SD, getTopX, topX=3), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
    rolledAssayData <- data.frame(rDT,row.names=1, check.names=F)
  }else if(method =="top1"){
    rDT <- DT[, lapply(.SD, getTopX, topX=1), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
    rolledAssayData <- data.frame(rDT,row.names=1, check.names=F)
  }

  # idScore
  if("idScore" %in% names(DT) ){

    ### get fData of highest scoring row per rollUP level
    indices <- DT[, .I[getMaxIndex(idScore) ], by=idx]$V1
    rolledFData <- data.frame(data.frame(DT)[indices,names(fData(eset))],row.names=rownames(rolledAssayData))
    # replace idScore colum, by summed score
    sumDT <- DT[, sum(idScore,na.rm=T), by=idx ]
    rolledFData$idScore = sumDT[,V1]
  }else{

    # @TODO inefficient solution -> speed up
    # allColumns but idx,intensity data and idScore
    #selColOld <- which(!(names(DT)  %in%  c(colnames(exprs(eset)),"idx","idScore")))
    selCol <- which((names(DT)  %in% names(fData(eset))))

    #rolledFDataOld <- data.frame(DT[, lapply(.SD, .getFirstEntry), by=idx, .SDcols=selCol],row.names=rownames(rolledAssayData))[,2:(length(selCol)+1)]
    rolledFData <- data.frame(DT[, lapply(.SD, .getFirstEntry), by=idx, .SDcols=selCol],row.names=rownames(rolledAssayData))[,names(fData(eset))]

    #		print(names(rolledFData))
    #		print("")
    #		print(names(rolledFDataOld))
  }

  ### concatenate allAccessions
  if("allAccessions" %in% names(DT)){
    rolledFData$allAccessions <- DT[, list( allAccessionsTMP  = paste(unique(unlist(strsplit(paste(allAccessions,collapse=";"),";"))),collapse=";") ), by = key(DT)]$allAccessionsTMP
  }

  rolledAssayData <- as.matrix(rolledAssayData)
  rolledAssayData[rolledAssayData == 0 ] <- NA
  #names(rolledFData) <- names(fData(eset))

  # reset isNormAnchor
  if(!is.null(rolledFData$isNormAnchor)){
    rolledFData$isNormAnchor <- T
  }
  # reset isFiltered
  if(!is.null(rolledFData$isFiltered)){
    rolledFData$isFiltered <- F
  }

  ### set peptides per protein
  rolledFData$nbPeptides <- getNbPeptidesPerProtein(eset)[as.character(rolledFData$proteinName)]

  ### if


  return(createExpressionDataset(expressionMatrix=rolledAssayData,expDesign=pData(eset),featureAnnotations=rolledFData))
}

#' Per Feature Normalization
#' @param eset ExpressionSet
#' @param normFactors matrix normalization factors (logged) (row names are proteins)
#' @return ExpressionSet object
#' @export
#' @details Example Usage: Normalize phospho peptide signals for Protein Changes
#' @note  No note
#' @references No references
#' @examples print("No examples")
perFeatureNormalization <- function(eset,normFactors){

  coveredPeptideSel <- fData(eset)$proteinName %in% rownames(normFactors)

  if(sum(coveredPeptideSel) == 0){
    warning("No shared entries between target data and normalisation data")
    return(eset)
  }

  # make sure target data and norm data have same dimensions
  if( !all(colnames(normFactors) %in% pData(eset)$condition) ){
    stop("Invalid norm factors")
  }

  # normalise log ratios
  exprs(eset)[coveredPeptideSel,]	<- exprs(eset)[coveredPeptideSel, ] - normFactors[as.character(fData(eset)[coveredPeptideSel,]$proteinName),pData(eset)$condition]


  return(eset)

}

#> 	pData(eset)
#condition isControl subject
#A_rep_1         A     FALSE       1
#A_rep_2         A     FALSE       2
#B_rep_1         B     FALSE       1
#B_rep_2         B     FALSE       2
#C_rep_1         C      TRUE       1
#C_rep_2         C      TRUE       2
#' Create Paired Expdesign
#' @param eset ExpressionSet
#' @return ExpressionSet object
#' @import Biobase
#' @export
#' @note  No note
#' @details  Add subject colum to phenoData design data.frame
#' @references NA
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
createPairedExpDesign  <-function(eset){

  # make sure all conditions include the same number of runs
  nbRunsPerCond <- unique(table(as.character(pData(eset)$condition)))

  if(length(nbRunsPerCond) != 1){
    print(pData(eset))
    stop("ERROR: createPairedExpDesign failed to create paired expDesign")
  }

  # Paired designnot meaningful if no replicates
  if(table(pData(eset)$condition)[1] == 1){
    cat("WARNING: createPairedExpDesign Paired experimental design not meaningful when no replicates\n" )
    return(eset)
  }

  pData(eset) <- cbind(pData(eset),subject=as.factor(rep(1:nbRunsPerCond,length(unique(pData(eset)$condition)))))

  return(eset)
}


#' Perform statistical test (mderated F-test)
#' @param eset ExpressionSet
#' @param adjust TRUE/FALSE adjust for multiple testing using Benjamini & Hochberg  (1995) method
#' @param log T/F log-transform expression values
#' @return list of pvalues
#' @export
#' @import limma Biobase
#' @importFrom stats model.matrix p.adjust
#' @note  No note
#' @details No details
#' @references Empirical Bayes method, Smyth (2004), \url{http://www.ncbi.nlm.nih.gov/pubmed/16646809}
#' @seealso \code{\link[limma]{eBayes}}
#' @examples print("No examples")
getFTestPValue = function(eset, adjust=F, log=T,...){

  pValues = list()

  # if no condition has replicates return NA's
  if(max(table(pData(eset)$condition)) == 1){
    pValues = rep(NA,nrow(eset))
    names(pValues) = rownames(eset)
    return(pValues)
  }

  if(log)(exprs(eset) <- log2(exprs(eset)))

  fit =eBayes(lmFit(eset,model.matrix(~condition, data=pData(eset)) ), 
              ...)[,-1]
  #pValues = list()
  pValues = fit$F.p.value
  names(pValues) = rownames(eset)

  if(adjust){ ### adjust for multiple testing using Benjamini & Hochberg (1995) method
    pValues <-p.adjust(pValues,method="BH")
  }

  return(pValues)
}


#' Roll up feature intensites per unique colum combination
#' @param eset ExpressionSet
#' @param featureDataColumnName vector of column names e.g. peptide or proteinName
#' @param method "sum", "mean" or "top3"
#' @return ExpressionSet object
#' @export
#' @details featureDataColumnName = c("peptide","charge","ptm"), method= c("sum"), sums up intensities per peptie modification charge state
#' @import Biobase
#' @note  No note
#' @references No references
#' @examples print("No examples")
rollUp <- function(eset, method = "sum", 	featureDataColumnName =  c("proteinName") ){

  ### apply filter
  eset <- eset[!fData(eset)$isFiltered,]

  # create rollup index (row names)
  if(length(featureDataColumnName) > 1 ){
    rnames = do.call(paste, c(fData(eset)[featureDataColumnName] , sep="_"))
  }else{
    rnames = fData(eset)[,featureDataColumnName]
  }


  # rollup feature data

  # add column of zero scores if no scores
  if(!("idScore" %in% names(fData(eset)))) fData(eset)$idScore = 0

  # add column of ones if nbFeatures doesn't exist
  if(!("nbFeatures" %in% names(fData(eset)))) fData(eset)$nbFeatures = 1

  # add index (integers)  and rnames columns to fData
  fData(eset) = cbind(fData(eset),index=1:nrow(eset),rnames)
  fdSummary = group_by( fData(eset), rnames) %>% summarise(idx=index[max(idScore) == idScore][1]
                                                           ,idScore = sum(idScore,na.rm=T)
                                                           ,allAccessions = paste(proteinName[!duplicated(proteinName)] ,collapse=";")
                                                           ,nbFeatures = sum(nbFeatures)
  )

  ### get fData of highest scoring row per rollUP level (do not include NA_IMP colummns)
  rolledFData = data.frame(fData(eset)[fdSummary$idx, !(colnames(fData(eset)) %>% grepl("^NA_IMP_[CI]NT\\.",.)) ],  row.names=fdSummary$rnames)
  # replace idScore column, by summed score, drop idScore column if all 0
  if(!all(fdSummary$idScore == 0)) rolledFData$idScore=fdSummary$idScore
  rolledFData$allAccessions=fdSummary$allAccessions
  rolledFData$nbFeatures=fdSummary$nbFeatures

  ### set peptides per protein
  #rolledFData$nbPeptides <- getNbPeptidesPerProtein(eset)[as.character(rolledFData$proteinName)]
  tmp = getNbPeptidesPerProtein(eset)
  rolledFData$nbPeptides <- tmp[match(as.character(rolledFData$proteinName),names(tmp))]
  
  # drop index and rnames columns
  rolledFData = rolledFData[, !(names(rolledFData) %in% c("index","rnames"))   ]

  # reset isNormAnchor
  if(!is.null(rolledFData$isNormAnchor)){
    rolledFData$isNormAnchor <- T
  }
  # reset isFiltered
  if(!is.null(rolledFData$isFiltered)){
    rolledFData$isFiltered <- F
  }

  # rollup expression data
  df = data.frame(exprs(eset),rnames )
  if(method =="sum"){

    # add NA intensties if tracked in fData
    df = data.frame(df,  getImputationMatrix(eset), getImputationMatrix(eset,method="count") )
    eRP = summarise_all(group_by_(df , .dots="rnames"  ) , funs(sum(.,na.rm=T))  )

    # add imputed intensities to rolledFData and rolled intensites to eset
    eRPCol = 1:(ncol(eset)+1)
    mImp = eRP[ !((1:ncol(eRP)) %in% eRPCol)]
    # set NA_IMP colnames
    colnames(mImp) = c(paste0("NA_IMP_INT.",colnames(eset)), paste0("NA_IMP_CNT.",colnames(eset)))

    #rolledFData = cbind(rolledFData,NA_IMP_INT=mImp[match(rownames(rolledFData),eRP$rnames), ]  )
    rolledFData = cbind(rolledFData,mImp[match(rownames(rolledFData),eRP$rnames), ]  )
    eRP = eRP[, eRPCol]

  }else if(method =="median"){
    eRP = summarise_all(group_by_(df , .dots="rnames"  ) , funs(median(.,na.rm=T))  )
  }else if(method =="mean"){
    eRP = summarise_all(group_by_(df , .dots="rnames"  ) , funs(mean(.,na.rm=T))  )
  }else if(method =="top3"){
    eRP = summarise_all(group_by_(df , .dots="rnames"  ) , funs(getTopX(.))  )
  }else if(method =="top1"){
    eRP = summarise_all(group_by_(df , .dots="rnames"  ) , funs(getTopX(.,topX = 1))  )
  }
  # store index as rownames
  eRP =  data.frame( eRP[-1], row.names=eRP[1] %>% unlist )
  colnames(eRP) = colnames(eset)

  # "X005_MLN.3.DDA_C17" -> "005_MLN-3-DDA_C17" issue
  names(eRP) = rownames(pData(eset))
  return(createExpressionDataset(expressionMatrix=eRP %>% as.matrix ,expDesign=pData(eset),featureAnnotations=rolledFData[match(rownames(eRP),rownames(rolledFData)),] ))

}


#' Impute missing values
#' @param eset ExpressionSet
#' @param method c("knn","ppca","gMin",lMin)
#' @param rowmax The maximum percent missing data allowed in any row to apply ppca and knn (if more missing values impute gmin). default 0.3
#' @return ExpressionSet
#' @export
#' @import impute
#' @importFrom  pcaMethods pca
#' @note  No note
#' @details
#' \itemize{
#' \item{gMin: half global minimum (0.1 percentile)}
#' \item{lMin: half local minimum}
#' \item{gMean: half global mean}
#' \item{lMean: half local mean}
#' \item{knn: Nearest neighbour averaging, as implemented in the impute::impute.knn function}
#' \item{ppca: An iterative method using a probabilistic model to handle missing values, as implemented in the pcaMethods::pca function.}
#' }
#' @references Accounting for the Multiple Natures of Missing Values in Label-Free Quantitative Proteomics Data Sets to Compare Imputation Strategies, Lazar et al (2016), \url{http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00981}
#' @seealso No note
#' @examples print("No examples")
sqImpute = function(eset,method="gmin", rowmax=0.3){

  esetImp = eset

  isNA = is.na(exprs(esetImp))
  if(sum(isNA) == 0) return(eset)

  if(ncol(esetImp) < 6){
    method="gmin"
    cat("INFO: sqImpute: method set to 'gmin'. nb. runs < 6\n")
  }

  # ignore case
  method = tolower(method)

  # some imputation method will fail if all entries per rowexp()

  #exprs(esetImp) = log2(exprs(esetImp))

  # get value at 0.5%
  set.seed(2018)
  gmin = (sample(exprs(esetImp),min(nrow(esetImp)*ncol(esetImp),10000))  %>% quantile(seq(0,1,0.001),na.rm=T))[6]

  if(method == "gmin"){

    exprs(esetImp)[is.na(exprs(esetImp))] = 0
    exprs(esetImp) = exprs(esetImp) + gmin

  }else if(method == "ppca"){
    suppressWarnings(suppressPackageStartupMessages(require(pcaMethods,warn.conflicts = F)))
    cat("INFO: Imputing missing values using ppca \n")
    # log transform to avoid neg imputed values
    #exprs(esetImp) =  exprs(esetImp) %>% log
    pc <- pcaMethods::pca(esetImp, nPcs=2, method="ppca")
    esetImp <- asExprSet(pc, esetImp)
    #exprs(esetImp) =  exprs(esetImp) %>% exp
    # neg values not accepted
    exprs(esetImp)[  exprs(esetImp) < 0 ] <- 0
    # add gmin
    exprs(esetImp) = exprs(esetImp) + gmin

  }else if(method == "knn"){
    suppressWarnings(suppressPackageStartupMessages(require(impute,warn.conflicts = F)))
    cat("INFO: Imputing missing values using knn \n")
    #exprs(esetImp) =  exprs(esetImp) %>% log
    invisible(capture.output(exprs(esetImp) <-  impute::impute.knn(exprs(esetImp), maxp=30000, rowmax=rowmax)$data))
    # add gmin
    exprs(esetImp) = exprs(esetImp) + gmin
    #exprs(esetImp) =  exprs(esetImp) %>% exp
  }else if(method == "lmin"){
    rowImp = apply(exprs(esetImp),1,min,  na.rm=T)/2
    exprs(esetImp) = exprs(esetImp) -  rowImp
    exprs(esetImp)[  is.na(exprs(esetImp)) ] <- 0
    exprs(esetImp) = exprs(esetImp) +  rowImp
  }else if(method == "gmean"){
    exprs(esetImp)[  is.na(exprs(esetImp)) ] <- mean(exprs(esetImp),  na.rm=T)
  }else if(method == "lmean"){
    rowImp = rowMeans(exprs(esetImp),  na.rm=T)
    exprs(esetImp) = exprs(esetImp) -  rowImp
    exprs(esetImp)[  is.na(exprs(esetImp)) ] <- 0
    exprs(esetImp) = exprs(esetImp) +  rowImp
  }else{
    stop("sqImpute: Unknown imputation method specified")
  }

  #exprs(esetImp) = 2^exprs(esetImp)

  # add gmin to rows with more than rowmax fraction NA's
   if(method %in% c("knn","ppca")){
       idxTooManyNA = (((is.na(exprs(eset)) %>% rowSums) / ncol(eset)) >= rowmax) %>% which
      for(i in idxTooManyNA){
        t = exprs(eset)[i,]
        t[is.na(t)] = 0
        exprs(esetImp)[i,] = t + gmin
        #t(rbind(t,exprs(esetImp)[i,])) %>% print
      }
      #cat("\t gmin method applied to ",signif(length(idxTooManyNA)/nrow(eset) * 100,3),"% of features\n"  )
    }

  # store imputed missing values
  impMatrix = exprs(esetImp)
  impMatrix[impMatrix != 0] = 0
  impMatrix[isNA] = exprs(esetImp)[isNA]
  fData(esetImp) = cbind(fData(esetImp), NA_IMP_INT=impMatrix, NA_IMP_CNT=isNA*1)

  return(esetImp)

}

#' Get matrix of imputed valeus in ExpressionSet matrix
#' @param eset ExpressionSet
#' @param method c("intensity","count") default intensity
#' @return matrix
#' @export
#' @note  No note
#' @details
#' @seealso No note
#' @examples print("No examples")
getImputationMatrix = function(eset, method="intensity"){

  # count or int
  regExpr = "^NA_IMP_INT\\."
  if(method == "count") regExpr = "^NA_IMP_CNT\\."

  isImpCol = grepl(regExpr,fData(eset) %>% colnames)

  if(sum(isImpCol) > 0){
    m = subset(fData(eset), select= isImpCol   )
    colnames(m) = colnames(eset)
  }else{
    m = exprs(eset)
    m[is.finite(exprs(eset)) | is.na(exprs(eset))] = 0
  }
  return(m)
}


#' Get fraction missing values per ratio/condition/run
#' @param eset ExpressionSet
#' @param method c("ratio","cond","run","count","intensity")
#' @return matrix
#' @export
#' @note  No note
#' @details
#' @seealso No note
#' @examples print("No examples")
getNAFraction = function(eset, method=c("ratio","intensity")){

  LEVELS = c("ratio","cond","run")
  SIGNALS = c("intensity","count")

  level = LEVELS[LEVELS %in% method]
  signal = SIGNALS[SIGNALS %in% method]

  if(signal == "count" ){
    exprs(eset) = getNbRolledUpFeatures(eset,method="matrix")
  }

  if(level == "run"){
    return((getImputationMatrix(eset, method = signal) / exprs(eset)) %>% as.matrix)
  }else{ # ratio or cond
    esetImp = eset
    exprs(esetImp) = getImputationMatrix(eset, method = signal) %>% as.matrix
    # medians
    intPerCond = getSignalPerCondition(eset,method ="sum" )
    naIntPerCond = getSignalPerCondition(esetImp, method="sum")

    if(level == "cond"){
      return((naIntPerCond/ intPerCond) %>% as.matrix)
    }else{ # ratio
      controlCond = eset$condition[eset$isControl][1]
      caseCond = eset$condition[!(eset$condition %in% controlCond)] %>% unique

      m = data.frame(row.names=rownames(eset))
      for(cond in caseCond){
        m = cbind(m,cond= (naIntPerCond[controlCond] + naIntPerCond[cond]) / (intPerCond[controlCond] + intPerCond[cond]))
      }
      names(m) = caseCond
      return(m %>% as.matrix )
    }
  }
}


#' Get number rolled up features per row
#' @param eset ExpressionSet
#' @param method c("vector","matrix") default vector
#' @return matrix
#' @export
#' @note  No note
#' @details
#' @seealso No note
#' @examples print("No examples")
getNbRolledUpFeatures = function(eset, method = "vector"){

  nbFeatures = rep(1,nrow(eset))
  if("nbFeatures" %in% colnames(fData(eset))  ){
    nbFeatures = fData(eset)$nbFeatures
  }

  if(method == "matrix"){
    nbFeatures = matrix(rep(nbFeatures,ncol(eset)),ncol=ncol(eset))
    colnames(nbFeatures) = colnames(eset)
    rownames(nbFeatures) = rownames(eset)
  }

  return(nbFeatures)
}


