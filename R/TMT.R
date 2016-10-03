# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

CALIBMIXRATIOS <- list()
CALIBMIXRATIOS [c(
				"sp|P00489|PYGM_RABIT"
				,"sp|P02789|TRFE_CHICK"
				,"sp|P01012|OVAL_CHICK"
				,"sp|P02666|CASB_BOVIN"
				,"sp|P00722|BGAL_ECOLI"
				,"tr|B6V3I5|B6V3I5_BOVIN"
		)] <- 1/c(4,0.25,2,0.5,2,0.5)

#CALIBMIXRATIOS <- subset(CALIBMIXRATIOS,!grepl("ECOLI",names(CALIBMIXRATIOS)))
#print(CALIBMIXRATIOS)

#' Get Thermo TMT impurity matrix
#' @param plexNb integer, 6 or 10 plex 
#' @return impurity matrix matrix
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getImpuritiesMatrix <- function (plexNb=6){
#getImpuritiesMatrix <- function (plexNb=6, test=F){
	
#	if(test){
#		
#		#### test from MSnbase
#		tmt126 <- c(0,0,6.1,0)
#		tmt127 <- c(0,0.5,6.7,0)
#		tmt128 <- c(0,1.1,4.2,0)
#		tmt129 <- c(0,1.7,4.1,0)
#		tmt130 <- c(0,1.6,2.1,0)
#		tmt131 <- c(0.2,3.2,2.8,0)
#		
#		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
#		rownames(reporterIonIsotopicDistrib) <- 1:6
#		
#	}else
	if(plexNb == 10){
		
		### FROM THERMO PRODUCT DATA SHEET
		
		# 10-plex
		tmt126 <- c(0,0,4.69,0)
		tmt127N <- c(0,0.4,6.5,0)
		tmt127C <- c(0,0.2,4.6,0.3)
		tmt128N <- c(0,0.9,4.7,0.2)
		tmt128C <- c(0.1,0.53,2.59,0)
		tmt129N <- c(0,0.73,2.49,0)
		tmt129C <- c(0,1.3,2.5,0)
		tmt130N <- c(0,1.2,2.8,2.7)
		tmt130C <- c(0.1,2.9,2.9,0)
		tmt131 <- c(0,2.36,1.43,0)
		
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0)
			,c(0,tmt127N,0),c(0,tmt127C,0)
			,c(0,tmt128N,0),c(0,tmt128C,0)
			,c(0,tmt129N,0),c(0,tmt129C,0)
			,c(0,tmt130N,0),c(0,tmt130C,0)
			,c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:10
		
	}else if(plexNb == 6){
		
		### FROM THERMO PRODUCT DATA SHEET
		# 6-plex old
#		tmt126 <- c(0,0,9.4,0.6)
#		tmt127 <- c(0.1,0.8,8.6,0.4)
#		tmt128 <- c(0.1,1.4,7.2,0.5)
#		tmt129 <- c(0.1,1.6,5.7,0.3)
#		tmt130 <- c(0.2,2.1,5.1,0.2)
#		tmt131 <- c(0.1,4.1,4.7,0.1)
		
#		tmt126 <- c(0,0,5.85,0.0)
#		tmt127 <- c(0.0,0.4,6.5,0)
#		tmt128 <- c(0.34,1.2,3.59,0)
#		tmt129 <- c(0,1.59,3.48,0)
#		tmt130 <- c(0.17,2.9,2.18,0)
#		tmt131 <- c(0.35,3.58,2.1,0)
		
		tmt126 <- c(0,   0  , 5.6  ,0)
		tmt127 <- c(0,   0.4, 5    ,0)
		tmt128 <- c(0,   1.2, 5.3  ,0.1)
		tmt129 <- c(0.2, 1.5, 4.1  ,0)
		tmt130 <- c(0.1, 2.6, 2.5  ,0)
		tmt131 <- c(0.1 ,3.3, 2.6  ,0)
		
		#### test from MSnbase
		#tmt126 <- c(0,0,6.1,0)
		#tmt127 <- c(0,0.5,6.7,0)
		#tmt128 <- c(0,1.1,4.2,0)
		#tmt129 <- c(0,1.7,4.1,0)
		#tmt130 <- c(0,1.6,2.1,0)
		#tmt131 <- c(0.2,3.2,2.8,0)
		
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
		
	}else if(plexNb == -6){	 ### test
		
		tmt126 <- c(0,   0  , 0  ,0)
		tmt127 <- c(0,   0, 0    ,0)
		tmt128 <- c(0,   0, 0  ,0)
		tmt129 <- c(0, 0, 4.1  ,0)
		tmt130 <- c(0, 0, 0  ,0)
		tmt131 <- c(0 ,3.3, 0  ,0)
		#tmt131 <- c(0 ,0, 0  ,0)
		
		plexNb <- 6
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
	
	}else if(plexNb == -60){	 ### test
		
		tmt126 <- c(0,   0, 0  ,0)
		tmt127 <- c(0,   0, 5  ,0)
		tmt128 <- c(0,   0, 0  ,0)
		tmt129 <- c(0, 1.5, 4.1,0)
		tmt130 <- c(0, 0, 0  ,0)
		tmt131 <- c(0 ,3.3, 0  ,0)
		#tmt131 <- c(0 ,0, 0  ,0)
		
		plexNb <- 6
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
		
	}
	
	###convert  reporterIonIsotopicDistrib -> impuritiesMatrix
	diagonals <- 1-apply(reporterIonIsotopicDistrib,1,sum)
	impuritiesMatrix <- diag(diagonals)
	
	for(i in 1:plexNb){
		
		### get affected channels
		affectedChannels <- c(-3,-2,-1,1,2,3) + i
		
		for(j in 1:6){
			affectedChannel <- 	affectedChannels[j]	
			if((affectedChannel > 0) & (affectedChannel <= plexNb)){
				impuritiesMatrix[i,affectedChannel] <- reporterIonIsotopicDistrib[i,j]
			}
		}
	}
	
	
	
	#diag(impuritiesMatrix) <- 1
	impuritiesMatrix <- t(impuritiesMatrix) 
	return(impuritiesMatrix)
	
}


#' Correct channel intensities based on Reporter ion Isotopic Distributions 
#' @param tmtData data.frame containing tmt channel intensities
#' @param impurityMatrix correction matrix
#' @return data.frame of corrected tmt intensities
#' @export
#' @note  No note
#' @details Same method as MSnbase, and described in Breitwieser et al. 2012 (Book Chapter)
#' @references NA 
#' @examples print("No examples")
purityCorrectTMT <- function(tmtData,impurityMatrix=impurityMatrix ){
	tmtDataCorrected <- matrix(nrow=nrow(tmtData),ncol=ncol(impurityMatrix))
	### solve linear system (spectrum by spectrum)
	for(i in 1:nrow(tmtData)){
		correctedSpectrumIntensities <- solve(impurityMatrix, tmtData[i,])
		tmtDataCorrected[i,1:ncol(impurityMatrix)] <- correctedSpectrumIntensities
	}
	
	tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- NA
	
#	if(invalidReplace  == "allZero" ){
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- 0
#	}else if(invalidReplace  == "allNA"){
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- NA
#	}else if(invalidReplace  == "allOrg"){ ### keep original values
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected <= 0)] <- tmtData[!is.na(tmtDataCorrected) & (tmtDataCorrected <= 0)] 
#	}### else Negative values are allowed
	
	colnames(tmtDataCorrected) <- colnames(tmtData)
	
	return(tmtDataCorrected)
}


#' Create Experimental Design
#' @param tag user input tag e.g. 1,2,3:4,5,6 indicating two condition with 3 reps each
#' @param nbPlex tmt 6  or 10 plex
#' @return expDesign data.frame
#' @export
#' @details The first listed condition is always the control condition
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
createExpDesign <- function(tag,nbPlex){
	
	sampleOrder <- as.numeric(unlist(strsplit(tag,"[\\,\\:]")))
	
	# make sure no duplicates, withing range etc.
	if(is.na(sampleOrder[1]) 
			| (max(table(sampleOrder))>1) 
			| (length(sampleOrder) > 10)
			| (max(sampleOrder) > 10)
			| (min(sampleOrder) < 1)
			| (length(sampleOrder) > nbPlex)
			){
		cat("ERROR: getExpDesign, INVALID EXPERIMENTAL DESIGN",tag,"\n")
		quit("no")
	}
	
	expDesign <- data.frame(row.names=sampleOrder,condition=rep("Ctrl",length(sampleOrder)),isControl=rep(F,length(sampleOrder) ))
	
	# create vector describing condition grouping e.g. c(3,3) ctrl 3 samples, 3 samples cond1
	condGrouping <- c()
	for(cond in unlist(strsplit(tag,":"))){
		condGrouping <- c(condGrouping,length(unlist(strsplit(cond,","))))
	}	
	
	expDesign$condition <- as.character(expDesign$condition)
	
	#browser()
	# update expDesign data.frame
	condNb <- 1
	sampleStartIdx <- 1 
	for(sampleCount in condGrouping ){
		
		if(condNb == 1){
			expDesign$condition[1:sampleCount] <- rep("Ctrl",sampleCount)
			expDesign$isControl[1:sampleCount] <- rep(T,sampleCount)
		}else{
			expDesign$condition[sampleStartIdx:(sampleStartIdx+sampleCount-1)] <- rep(paste("cond",condNb-1,sep="_"),sampleCount)
			expDesign$isControl[sampleStartIdx:(sampleStartIdx+sampleCount-1)] <- rep(F,sampleCount)
		}
		condNb <- condNb+1
		sampleStartIdx <-  sampleStartIdx+sampleCount
	}	
	
	return(expDesign)
}


#' Sum up raw intensities per protein and channel. keep track of number of summed spectra and unique peptides
#' @param intData data.frame of intensities per channel
#' @param proteinACs vector of protein accession numbers
#' @param peptides vector of peptide sequneces
#' @param minNbPeptPerProt minimal number of peptides per protein
#' @return list containing  3 objects 1) data.frame of channel intensities per protein ac, 2) vector listing number of summed spectra per protein, 3) vector listing number of summed peptides per protein
#' @export
#' @details NA
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getIntSumPerProtein  <- function(intData,proteinACs,peptides,minNbPeptPerProt=1){
	uniqueProteinACs <- unique(proteinACs)
	nbProteins <- length(uniqueProteinACs)
	
	perProteinIntSum <- matrix(nrow=nbProteins, ncol=ncol(intData))
	
	c <-1
	### init progressio bar
	pbSum <- txtProgressBar(min = 0, max = nbProteins, style = 3)
	
	spectraPerProtein <- c()
	peptidesPerProtein <- c()
	for(ac in uniqueProteinACs){
		
		# disp/increment progression bar
		setTxtProgressBar(pbSum, c)
		
		
		sel <- proteinACs %in% ac 
		sset 	<- intData[sel,]
		uPept <- unique(peptides[sel])
		
		if(sum(sel) > 1){
			perProteinIntSum[c,] <- apply(sset,2,sum)
		}else{
			perProteinIntSum[c,] <- as.vector(unlist(sset))
		}
		
		spectraPerProtein <- c(spectraPerProtein,nrow(sset))
		peptidesPerProtein <- c(peptidesPerProtein,length(unique(peptides[sel])))
		
		c <- c+1
	}
	# terminate progession bar
	close(pbSum)
	cat("\n")
	
	perProteinIntSum  <- data.frame(perProteinIntSum)
	rownames(perProteinIntSum) <- uniqueProteinACs
	names(perProteinIntSum) <- paste("channel_",1:ncol(intData),sep="")
	
	### FILTER based on min peptides per protein
	ok <- peptidesPerProtein >= minNbPeptPerProt
	if(sum(ok) < 2 ){
		cat("ERROR: getIntSumPerProtein. Not enough Features, Review filtering parameter minNbPeptPerProt ", minNbPeptPerProt,"\n")
		quit(status=-1)
	}
		
	perProteinIntSum <- perProteinIntSum[ok,]
	spectraPerProtein <- spectraPerProtein[ok]
	peptidesPerProtein <- peptidesPerProtein[ok]

	ret <- list()
	ret$perProteinIntSum <- perProteinIntSum
	ret$spectraPerProtein <- spectraPerProtein
	ret$peptidesPerProtein <- peptidesPerProtein
	
	return(ret)
	
}


#' @export
.intensityAdjustment <- function(eset,esetCalibMix){
	
	#
	noiseFractionIQRThrs <- 0.4
	#samplePairVariationThrs <- 1
	minAdjustedIntThrs <- 0.1
	
	###################### GET GLOBAL NOISE FRACTION  ######################
	
	# channel Intensity Fraction , non cal mix
	F_I_channel <- apply(exprs(eset),2,sum,na.rm=T) / sum(exprs(eset),na.rm=T)
	
	# intenstiy reporter ion cluster cal mix, I_CM (@TODO ignore NAs)
	I_CM <- apply(exprs(esetCalibMix),1,sum) ### THIS LEADS TO MORE NAs IN F_N_EST
	
	# I - reporter ion cluster non cal mix. (How to handle NAs ?)
	I <- apply(exprs(eset),1,sum,na.rm=T)
	
	F_N_EST <- data.frame(row.names=row.names(esetCalibMix))
	
	refRatio <- as.vector(unlist(CALIBMIXRATIOS[as.character(fData(esetCalibMix)$proteinName)]))
	
	for(i in c(1,3,5,7,9)){
		
		# variaiton in sanmple amounts of pais
#		tmp <- F_I_channel[c(i,i+1)]	
#		cv <- sd(tmp) / mean(tmp)
#		
#		# apply different cut-off for sample pairs e.g. [0.3,0.4,0.2,0.3,0.3] (per "i")
#		if(cv > samplePairVariationThrs){
#			F_N_EST <- cbind(F_N_EST, rep(NA,nrow(F_N_EST)))
#			
#			cat("WARN: IGNORED CAL-MIX PAIR ", i,":",i+1, " C.V. ", round(cv,2)*100,"% \n")
#			
#		}else{
		# 1) Calculate interfenece intensity per channel pair A (N_CM_Pair)
		# N_CM_Pair =  (I_CM_Pair_1 - R_Ref*I_CM_Pair_2)/ 	(1-R_Ref)	
		N_CM_Pair <- (exprs(esetCalibMix)[,i+1] - (refRatio*exprs(esetCalibMix)[,i]) )  /(1-refRatio)
		# set negative interference levels to zero
		N_CM_Pair[N_CM_Pair < 0] <- NA
		
		# 2) Calculate total intenference intensity in spectrum
		#	 N_CM = N_CM_Pair / (I_Pair / I_Tot)
		N_CM = (N_CM_Pair) / mean(F_I_channel[c(i,i+1)]) 
		
		# 3) interference fraction per reporter ion cluster
		#	 F_N = N_CM / I_CM
		F_N = N_CM / I_CM
		F_N[F_N > 1] <- NA
		
		F_N_EST <- cbind(F_N_EST,F_N)
#		}
	}
	names(F_N_EST) <- 1:ncol(F_N_EST)
	
	#F_N_GLOBAL <- median(apply(F_N_EST,2,median,na.rm=T),na.rm=T)
	selectedPairs <- apply(F_N_EST,2,IQR,na.rm=T) < noiseFractionIQRThrs
	F_N_GLOBAL <- min(apply(F_N_EST,2,median,na.rm=T)[selectedPairs],na.rm=T)
	
	
	if(is.na(F_N_GLOBAL)){
		stop("Large Differnces in sample statrting amounts. Ratio Adjustment not possible")
	}
	
	###################### GET GLOBAL NOISE FRACTION END ######################
	
	# Calculate Total interference intensity per reporter ion cluster (Non-cal mix data)
	# N = I * F_N_GLOBAL
	N = I * F_N_GLOBAL
	
	# Calculate channel  Interference intensity according to Channel Intensity Fractions (Non-cal mix data)
	# N_channel = N * F_I_channel
	N_channel = N * matrix(rep(F_I_channel,each=nrow(eset)),nrow=nrow(eset))
	
	### temp set NA's to Inf (to allow for implementation of minAdjustedIntThrs)
	I_channel <-  exprs(eset)
	I_channel[is.na(I_channel)] <- Inf
	
	# calcualte adjusted channel intensity
	# I_Adj_channel = I_channel	- N_channel 
	# if  I_Adj_channel < 0.05* I_channel
	#	Then  I_Adj_channel = 0.05* I_channel
	I_Adj_channel = I_channel - N_channel 
	I_Adj_channel[I_Adj_channel < (minAdjustedIntThrs* I_channel)] <-  minAdjustedIntThrs* I_channel[I_Adj_channel < (minAdjustedIntThrs * I_channel)]
	I_Adj_channel[I_Adj_channel == Inf] <- NA
	
	esetAdj <- eset
	exprs(esetAdj) <- I_Adj_channel
	
	########################## ADJUST CAL MIX #############################
	I_CM <- apply(exprs(esetCalibMix),1,sum,na.rm=T)
	N_CM <- I_CM * F_N_GLOBAL
	N_CM_channel = N_CM * matrix(rep(F_I_channel,each=nrow(esetCalibMix)),nrow=nrow(esetCalibMix))
	
	I_CM_channel <-  exprs(esetCalibMix)
	I_CM_channel[is.na(I_CM_channel)] <- Inf
	
	I_CM_Adj_channel = I_CM_channel - N_CM_channel 
	I_CM_Adj_channel[I_CM_Adj_channel < (minAdjustedIntThrs* I_CM_channel)] <-  minAdjustedIntThrs* I_CM_channel[I_CM_Adj_channel < (minAdjustedIntThrs * I_CM_channel)]
	I_CM_Adj_channel[I_CM_Adj_channel == Inf] <- NA
	
	esetCalMixAdj <- esetCalibMix
	exprs(esetCalMixAdj) <- I_CM_Adj_channel
	
	########################## ADJUST CAL MIX END ##########################
	
	ret <- list()
	ret$esetAdj <- esetAdj
	ret$esetCalMixAdj <- esetCalMixAdj
	ret$noiseFraction <- F_N_EST
	ret$globalNoiseFraction <- F_N_GLOBAL
	ret$selectedPairs <- selectedPairs
	
	return(ret)
}




### filter for calMix protein and add reference ratios to feature data
#' @export
.getCalibMixEset <- function(eset, nbPlex = 10){
	
	# select for claib mix proteins
	esetCalibMix <- eset[(fData(eset)$proteinName %in% names(CALIBMIXRATIOS)) ,]
	# add refernce ratios to feature data
	fData(esetCalibMix)$refRatio <- as.vector(unlist(CALIBMIXRATIOS[as.character(fData(esetCalibMix)$proteinName)]))
	
	if(nbPlex == 10){
		# 1:2 20	5:6 4
		# 3:4 100	1:2,7:8  20
		# 5:6 4	--> 3:4,9:10 100
		# 7:8 20	
		# 9:10 100
		# 
		pData(esetCalibMix)$condition <- as.factor(paste("cond",c(1,2,3,4,5,6,1,2,3,4),sep="_"))
		pData(esetCalibMix)$isControl <- rep(F,10)
		pData(esetCalibMix)$isControl <- c(T,F,F,F,F,F,T,F,F,F)
	}else{
		stop("ERROR: Unkwon nbPLex option", nbPlex)
	}
	
	return(esetCalibMix)
}

### filter for calMix protein and add reference ratios to feature data
#' @export
.getCalibMixPairedEset <- function(esetCalibMix){
	
	# make sure input eset only contain claibration mix proteins
	esetCalibMix <- esetCalibMix[fData(esetCalibMix)$proteinName %in% names(CALIBMIXRATIOS),]
	
	# pair up featureData
	esetPair <- esetCalibMix[,1:2]
	
	# add refernce ratios to feature data in case it's not already there
	if(!("refRatio" %in% names(fData(esetCalibMix)))){
		fData(esetCalibMix)$refRatio <- as.vector(unlist(CALIBMIXRATIOS[as.character(fData(esetCalibMix)$proteinName)]))
	}
	
	fData(esetPair) <- rbind(fData(esetPair) 
			, fData(esetPair)
			, fData(esetPair)
			, fData(esetPair)
			, fData(esetPair)
	)
	
	# add calibration mix dilution tag to featureData
#	calMixDilution <- c(rep(20,nrow(esetPair))
#			,rep(100,nrow(esetPair))
#			,rep(4,nrow(esetPair))
#			,rep(20,nrow(esetPair))
#			,rep(100,nrow(esetPair))
#	)
	calMixDilution <- c(rep(20,nrow(esetPair))
			,rep(100,nrow(esetPair))
			,rep(4,nrow(esetPair))
			,rep(20,nrow(esetPair))
			,rep(4,nrow(esetPair))
	)
	
	fData(esetPair) <- cbind(fData(esetPair),calMixDilution=calMixDilution)
	
	#add dilution to peptide and proteinName
	fData(esetPair)$proteinName <- paste(fData(esetPair)$proteinName,"_",calMixDilution,sep="") 
	fData(esetPair)$peptide <- paste(fData(esetPair)$peptide,"_",calMixDilution,sep="") 
	
	# pair up expression data
	exprs(esetPair) <- rbind(
			exprs(esetPair)[,1:2]
			,exprs(esetCalibMix[,3:4])
			,exprs(esetCalibMix[,5:6])
			,exprs(esetCalibMix[,7:8])
			,exprs(esetCalibMix[,9:10])
	)
	rownames(esetPair) <- 1:nrow(esetPair)
	
	# get rid of globalNormFactors column
	pData(esetPair) <- pData(esetPair)[,1:2]
	
	return(esetPair)
	
}

##' Get linear model explaning log2 ratio as a function of log2 tmt ratio 
##' @param eset paired calibration mix
##' @return linear model
##' @export
##' @note  No note
##' @details Uses linear model of log tmt ratio vs  log ref ratio
##' @references NA 
##' @examples print("No examples")
#getRatioCorrectionFactorModel <- function(eset){
#	log2RefRatio <- log2(fData(eset)$refRatio)
#	log2TmtRatio <- getRatios(eset)[,1]
#	ok <- is.finite(log2RefRatio) & is.finite(log2TmtRatio)
#	data <- data.frame(refRatio=log2RefRatio[ok],tmtRatio=log2TmtRatio[ok])
#	
#
#	
#	fit <- lm(  refRatio ~ tmtRatio  ,data=data)
#	if(abs(fit$coefficients[1]) > 0.1){
#		cat("WARN: getRatioCorrectionFactorModel: Large intercept\n")
#		#print(fit)
#	}
#	
#	#fit$coefficients[1] <- 0
#	fit <- lm( refRatio ~ 0 + tmtRatio  ,data=data)
#	
#	return(fit)
#	
#}

# Algorithm
# N_CM_Pair - Noise in calmix pair 
# N_CM - Tot noise in cal mix tmt cluster
# F_N - Fraction Noise in tmt cluster
# F_N_GLOBAL - Global estimate of Fraction Noise in tmt cluster
# I_Pair 
# I_CM_Pair - intenstiy cal mix pair
# I_CM - intenstiy tot cal mix cluster
# R_Ref - refernce ratio
# I - intensity total   
# N - total interfernece intensity
# F_I_channel - channel fraction intensity
# N_channel - channel noise
# I_channel	- channel intensity
# I_Adj_channel - adjusted channel intensity

# Each channel pair of all cal mix spectra gives an estimate of the globl cluster interference fraction 
# Large diffenerces in channel intensity in a channel pair is problemeatic
# Only use channel pairs where sd(I_channel) / mean(I_Pair) < 15% 

# Estimate interference fraction per reporter ion cluster (using calmix data and non-cal mix data)
# 1) Calculate interfenece intensity per channel pair A (N_CM_Pair)
# N_CM_Pair =  (I_CM_Pair_1 - R_Ref*I_CM_Pair_2)/ 	(1-R_Ref)	
# 2) Calculate total intenference intensity in spectrum	
# N_CM = N_CM_Pair / (I_Pair / I_Tot)
# 3) interference fraction per reporter ion cluster
# F_N = N_CM / I_CM

# Calculate Total interference intensity per reporter ion cluster (Non-cal mix data)
# N = I *  F_N_GLOBAL
# Calculate channel  Interference intensity according to Channel Intensity Fractions (Non-cal mix data)
# N_channel = N * F_I_channel
# chalcualte adjusted channel intensity
# I_Adj_channel = I_channel	- N_channel 
# if  I_Adj_channel < 0.05* I_channel
#	Then  I_Adj_channel = 0.05* I_channel


##' Adjust TMT intensitis discarding estimated intenference (Experimental)
##' @param eset 10 plex TMT dataset
##' @param eset calibration mix
##' @return intensityAdjustment list 
##' @note  No note
##' @details  Adjust TMT intensitis discarding estimated intenference (Experimental)
##' @references NA 
##' @examples print("No examples")



