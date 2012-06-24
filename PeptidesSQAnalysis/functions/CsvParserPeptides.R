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
	
	### fix expDesignDF
	
	expDesignDF <- expDesignDF[, names(expDesignDF) != "Bestpeptidematch"]
	expDesignDF[,length(names(expDesignDF))] <- expDesignDF[,length(names(expDesignDF))] -1
	
	### get number of samples
	nbSamples <- 0
	for(cond in names(expDesignDF)){
		nbSamples <- nbSamples+ expDesignDF[,cond]
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
	
	#expDesignObj$expDesign <- expDesign
	expDesignObj$expDesign <- expDesignDF
	expDesignObj$nbSamples <- nbSamples
	#expDesignObj$condition.names <- condition.names
	#expDesignObj$orgCondNames <- orgCondNames
	expDesignObj$isEnoughReplicates <- isEnoughReplicates
	
	return(expDesignObj)
	
}


filterProgOutDataPeptides <- function(progOutData, verbose=verbose, selectedProteinName=selectedProteinName, selectedModifName=selectedModifName,precursorMassFilter=precursorMassFilter, fdrCutoff=0.01  ){
	
	if(verbose){
		print("FILTERING PROTEIN LIST")
	}
	
	### FILTER PROTEIN LIST (NAME,PEPTIDE COUNTS)
	
	### filter by proteins name
	proteinNameFilterCond <- regexpr(selectedProteinName,progOutData$Accession,ignore.case=TRUE) > -1
	
	### remove contaminants
	nonConProteinFilter  <- !isCon(progOutData$Accession)
	#peptidesForQuantFilterCond <- progOutData$Peptides.used.for.quantitation >= minNbPeptidesPerProt
	
	### filter by modification name 	
	modifNameFilterCond <- rep(1,length(progOutData$modifs))
	if(selectedModifName != "." ){
		modifNameFilterCond <- regexpr(selectedModifName,progOutData$modif,ignore.case=TRUE) > -1
	}
	proteinsBeforeFilter <- nrow(progOutData)
	progOutData <- progOutData[proteinNameFilterCond & nonConProteinFilter & modifNameFilterCond,]
	proteinsAfterFilter <- nrow(progOutData)
	
	if(verbose){
		print(paste("    ", proteinsBeforeFilter- proteinsAfterFilter," peptides filtered out"))
	}
	
	
	### FILTER PROTEIN LIST (NAME,PEPTIDE COUNTS) END
	
	### precursor mass filter
	minMassDiffs <- progOutData$minMassDiffs
	
	### for precursor mass plot
	massDiffPlotObj <- list()
	massDiffPlotObj$minMassDiffs <- minMassDiffs
	massDiffPlotObj$scores <- progOutData$Confidence.score
	massDiffPlotObj$isDecoy <- isDecoy(progOutData$Accession)
	
	peptidesBeforeFilter <- nrow(progOutData)
	progOutData <- progOutData[!(minMassDiffs > precursorMassFilter),]
	
	if(verbose){
		peptidesAfterFilter <- nrow(progOutData)
		print("precursorMassFilter")
		print(paste("    ", peptidesBeforeFilter - peptidesAfterFilter," peptides filtered out"))
	}
	
	### CALCULATE FDR

	
	###FOR PLOTS
	decoyCond <- isDecoy(progOutData$Accession)
	scores <- progOutData$Confidence.score
	if(verbose){
		print("---- FDR CALCULATION ----")
	}
	qvals <- getIdLevelQvals(scores,decoyCond)
	
	### CALCULATE FDR END
	
	### FILTER PROTEIN LIST (FDR, DECOY)
	
	if(verbose){
		print("---- FILTER VALID PROTEIN IDS ----")
	}
	
	fdrFilterCond <- qvals < fdrCutoff

	proteinsBeforeFilter <- nrow(progOutData)
	progOutData <- progOutData[!decoyCond & fdrFilterCond,]
	
	if(verbose){
		proteinsAfterFilter <- nrow(progOutData)
		print("Id. level fdr filter")
		print(paste(proteinsBeforeFilter,"    ", proteinsBeforeFilter- proteinsAfterFilter," peptides filtered out"))
	}
	
	### FILTER PROTEIN LIST (FDR, DECOY) END
	
	filterObj = list()
	
	filterObj$progOutData <- progOutData
	
	### for id rate related plots
	filterObj$massDiffPlotObj <- massDiffPlotObj
	filterObj$qvals <- qvals
	filterObj$decoyCond <- decoyCond
	filterObj$scores <- scores
	
	return(filterObj)
	
}










### sum up data over unique peptides -> CREATE DATA FRAME OF UNIQUE PEPTIDES (max score & min mass Error (ppm) selected )
getSummedUpPeptideProgOutData <- function(data,intStartColumn =intStartColumn,nbSamples=nbSamples, proteaseTarget=proteaseTarget,verbose=F  ){
	
	### discard non-peptide-annotated features 
	data <- data[nchar(as.character(data$Sequence)) > 0,] 
	
	### dicard Mass column (not present in all .csv files, version issue?)
	data <- data[,!(names(data) == "Mass")] 
	
	rawIntColumnStart <- intStartColumn+nbSamples
	rawIntColumnEnd <- rawIntColumnStart+nbSamples-1
	spectralCountColumnsStart <-  rawIntColumnEnd +(2*nbSamples)+ 9	
	spectralCountColumnsEnd <- spectralCountColumnsStart+nbSamples-1
	
	petideIds <-  paste(data$Sequence,data$Variable.modifications...position..description. ,sep="::")
	uniquePeptideIds <- unique(petideIds)
	
	peptides <- c()
	modifs <- c()
	Accession <- c()
	Description =  c()
	minMassDiffs <- c()
	Confidence.score <-c()
	spc <- c()
	rawIntensities <- c()
	nbSummedEntries <- c()
	
	if(verbose){
		print("Grouping features by unique peptide sequence.. this can also take a while")
	}
	
	pbSum <- txtProgressBar(min = 0, max = length(uniquePeptideIds), style = 3)
	
	i <-1
	for(pept in uniquePeptideIds){
		
		setTxtProgressBar(pbSum, i)
		i <-  i+1
		
		indices = which(pept == petideIds )
		index1 = indices[1]
		
		peptides = c(peptides, as.character(data$Sequence[index1]))
		modifs = c(modifs,as.character(data$Variable.modifications...position..description.[index1]))
		Accession = c(Accession,as.character(data$Protein[index1]))
		Description =  c(Description,as.character(data$Description[index1]))
		minMassDiffs = c(minMassDiffs, min(data$Mass.error..ppm.[indices]))
		Confidence.score = c(Confidence.score,max(data$Score[indices]))
		nbSummedEntries = c(nbSummedEntries, length(indices))
		
		###sum up intensities
		intensities.tmp = data[indices,rawIntColumnStart:rawIntColumnEnd]
		sumRawIntensities = list()
		sumRawIntensities[[pept]] = apply(intensities.tmp,2,FUN=sum)
		rawIntensities = c(rawIntensities,sumRawIntensities) 
		
		###sum up spectral counts
		spC.tmp = data[indices,spectralCountColumnsStart:spectralCountColumnsEnd]
		sumSpC = list() 
		sumSpC[[pept]] = apply(spC.tmp,2,FUN=sum)
		spc = c(spc,sumSpC)
		
	}
	close(pbSum)
	cat("\n")
	
	### COUNT MIS-CLEAVAGE SITES
	proteaseRegexpr <- paste("[",proteaseTarget,"]",".",sep="")
	missedCleavageSites <- unlist(lapply(peptides,FUN= function(f){return(sum(unlist(gregexpr(proteaseRegexpr,f))>0))}))
	### COUNT MIS-CLEAVAGE SITES END
	
	rawData.tmp <- t(data.frame(rawIntensities))
	spectralCounts.tmp <- t(data.frame(spc)) 
	
	dataSummed <- data.frame(peptides
			,modifs
			,Accession
			,Description		
			,minMassDiffs
			,Confidence.score
			,missedCleavageSites		
			,nbSummedEntries
			,rawData.tmp
			,spectralCounts.tmp
			,row.names = uniquePeptideIds
	)
	
	
	
	return(dataSummed)
	
}
