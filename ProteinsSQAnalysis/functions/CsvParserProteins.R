# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

# A) extract experimental design
# B) x			extract data
# C) x 			filter raw data 
# D) create normdata object
# E) create spectral counts object

parseExpDesign <- function(progenesisFile, verbose=verbose){
	
	line2header <- read.csv(file=progenesisFile,head=TRUE,sep=",",check.names=FALSE, skip=1,nrows=1, strip.white=TRUE)
	line2header <- gsub(" ","",names(line2header))
	line2header.names.unique <- unique(line2header)
	condition.names <- line2header.names.unique[nchar(line2header.names.unique) > 0]
	
	for(cond in condition.names){
		line2header[line2header == cond] <- as.character(cond)
	}

	
	expDesign <- list()
	for(cond in condition.names){
		expDesign[[cond]] <- 0
	}
	
	currentLab <- names(expDesign)[1] 
	replicates <- 0
	
	for(lab in line2header){
		
		if(lab %in% names(expDesign)){
			expDesign[[currentLab]] <- replicates + 1
			currentLab <- lab
			replicates <- 0
		}else{
			replicates <- replicates+1
		}
	}
	
	### test, list 2 df
	expDesignDF <- data.frame(row.names="reps")
	for(i in names(expDesign)){
		expDesignDF <- data.frame(expDesignDF,expDesign[i])
	}
	names(expDesignDF) <- names(expDesign)
	
	### get number of samples
	nbSamples <- 0
	for(cond in condition.names){
		nbSamples <- nbSamples+ expDesign[[cond]]
	}
	
	isEnoughReplicates <- min(expDesignDF) > 1
	
	if(userOptions$verbose){
		print("Experimental design")
		print(expDesignDF)
	}	
	
	### READ EXPERIMENTAL DESIGN END
	
	if(userOptions$verbose){
		print("condition.names")
	}
	
	expDesignObj <- list()
	
	expDesignObj$expDesign <- expDesignDF
	expDesignObj$nbSamples <- nbSamples
	expDesignObj$isEnoughReplicates <- isEnoughReplicates
	
	return(expDesignObj)
	
}


filterProgOutDataProteins <- function(progOutData, verbose=verbose, selectedProteinName=selectedProteinName, minNbPeptidesPerProt=minNbPeptidesPerProt,fdrCutoff=0.01 ){
	
	if(verbose){
		print("FILTERING PROTEIN LIST")
	}
	
	### FILTER PROTEIN LIST (NAME,PEPTIDE COUNTS)
	
	proteinNameFilterCond <- regexpr(selectedProteinName,progOutData$Accession,ignore.case=TRUE) > -1
	nonConProteinFilter  <- !isCon(progOutData$Accession)
	peptidesForQuantFilterCond <- progOutData$Peptides.used.for.quantitation >= minNbPeptidesPerProt
	
	proteinsBeforeFilter <- nrow(progOutData)
	progOutData <- progOutData[proteinNameFilterCond & nonConProteinFilter & peptidesForQuantFilterCond,]
	proteinsAfterFilter <- nrow(progOutData)
	if(verbose){
		print(paste("    ", proteinsBeforeFilter- proteinsAfterFilter," proteins filtered out"))
	}
	
	
	### FILTER PROTEIN LIST (NAME,PEPTIDE COUNTS) END
	
	### CALCULATE FDR
	
	if(verbose){
		print("---- FDR CALCULATION ----")
	}
	
	### ALL
	decoyCond <- isDecoy(progOutData$Accession)
	scores <- progOutData$Confidence.score
	decoyScores <- scores[decoyCond]
	targetScores <- scores[!decoyCond]
	qvals <- getIdLevelQvals(scores,decoyCond)
	
	### CALCULATE FDR END
	
	### FILTER PROTEIN LIST (FDR, DECOY)
	
	if(verbose){
		print("---- FILTER VALID PROTEIN IDS ----")
	}
	
	fdrFilterCond <- qvals < fdrCutoff
	
	proteinsBeforeFilter <- nrow(progOutData)
	progOutData <- progOutData[!decoyCond & fdrFilterCond,]
	
	row.names(progOutData) <- progOutData$Accession
	
	proteinsAfterFilter <- nrow(progOutData)
	
	if(verbose){
		print(paste("    ", proteinsBeforeFilter- proteinsAfterFilter," proteins filtered out"))
	}
	
	### FILTER PROTEIN LIST (FDR, DECOY) END
	
	filterObj = list()
	
	filterObj$progOutData <- progOutData
	
	### for id rate related plots
	filterObj$qvals <- qvals
	filterObj$decoyCond <- decoyCond
	filterObj$scores <- scores
	
	return(filterObj)
	
}

