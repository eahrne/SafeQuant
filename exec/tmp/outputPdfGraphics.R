# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


######################################## GRAPHICS
cat("CREATING OUTPUT FILES\n")
pdf(userOptions$pdfFilePath)

if(userOptions$verbose){
	cat("CREATING PLOTS\n")
}

### COLORS
CONDITIONCOLORS <-  data.frame(COLORS[(1:EXPDESIGNOBJ$nbConditions)], row.names=names(EXPDESIGNOBJ$expDesign ))
names(CONDITIONCOLORS) <- c("cond")

#SAMPLECOLORS <-  data.frame(1:length(names(UNGROUPEDRAWINTDATA)), row.names=names(UNGROUPEDRAWINTDATA))
### COLORS END

### EXP DESIGN PLOT
if(userOptions$isDispExpDesign){
	plotExpDesign(EXPDESIGNOBJ$expDesign ,sampleNames=names(UNGROUPEDRAWINTDATA), controlCondition=CONTROLCONDITION,condColors=CONDITIONCOLORS, version=VERSION )
}
### EXP DESIGN PLOT END

### FDR PLOTS
if(exists("LFQDIR") && userOptions$isFdrPlots){
	
	if(ISPEPTIDEANALYSIS){
		
		par(mfrow=c(2,1))	
		#plotTargetDecoyRationPerMassDiffBin(FILTEROBJ$massDiffPlotObj$minMassDiff,FILTEROBJ$massDiffPlotObj$isDecoy)
		plotMassDiffDistibution(FILTEROBJ$massDiffPlotObj$minMassDiff,FILTEROBJ$massDiffPlotObj$isDecoy,precursorMassFilter=userOptions$precursorMassFilter)
		par(mar=c(5.1,4.1,0.1,4.1))		
		### plot precursor mass diff vs. score
		plotMassDiffVsScore(FILTEROBJ$massDiffPlotObj$minMassDiff,FILTEROBJ$massDiffPlotObj$scores,FILTEROBJ$massDiffPlotObj$isDecoy, precursorMassFilter=userOptions$precursorMassFilter)
		par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
		
		
		### id related plots
		par(mfrow=c(3,1))
		plotScoreDistrib(FILTEROBJ$scores[!FILTEROBJ$decoyCond],FILTEROBJ$scores[FILTEROBJ$decoyCond],nbBins=100,scoreName="Confidence score",ylab="Feature Count",title="")
		plotROC(FILTEROBJ$qvals[!FILTEROBJ$decoyCond],fdrMax=0.05,nbDataPoints=100,ylab="Valid Features",fdrCutOff=userOptions$fdrCutoff
				,title=paste("# Features: ",sum(FILTEROBJ$qvals < userOptions$fdrCutoff),", # Unique Peptides: ",nrow(PROGOUTDATA)," (FDR ",userOptions$fdrCutoff,")",sep="" )) 
		plotIdScoreVsFDR(FILTEROBJ$scores,FILTEROBJ$qvals, xlim=c(0,max(FILTEROBJ$scores)), xlab="Confidence score", ylab="FDR", type="l")
		par(mfrow=c(1,1))
		
		
	}else{
		### id related plots
		par(mfrow=c(3,1))
		plotScoreDistrib(FILTEROBJ$scores[!FILTEROBJ$decoyCond],FILTEROBJ$scores[FILTEROBJ$decoyCond],nbBins=100,scoreName="Confidence score",ylab="Protein Count",title="")
		plotROC(FILTEROBJ$qvals[!FILTEROBJ$decoyCond],fdrMax=0.05,nbDataPoints=100,ylab="Valid Proteins",fdrCutOff=userOptions$fdrCutoff) 
		plotIdScoreVsFDR(FILTEROBJ$scores,FILTEROBJ$qvals, xlim=c(0,max(FILTEROBJ$scores)), xlab="Confidence score", ylab="FDR", type="l")
		par(mfrow=c(1,1))
	}
	
	
}
### FDR PLOTS END

### INTENSITY DISTRIBUTION PLOTS
if(userOptions$isIntensityDistributionPlots){
	
	### plots all against all per condition
	for(cond in names(EXPDESIGNOBJ$expDesign )){
		pairs.annot(log(GROUPEDNORMINTDATA[[cond]]), main=cond)
	}
	
	### plot all against all condition median intensities
	pairs.annot(log(MEDIANNORMINTDATA))
	
	### two plots per slide
	par(mfrow=c(2,1))
	
	### boxplot of CVS per condition
	if(exists("CVSPERCOND")){ ### does not exists unless replicates
		boxplot(data.frame(CVSPERCOND[names(EXPDESIGNOBJ$expDesign )])*100,ylab= "C.V. (%)" ,col=as.character(CONDITIONCOLORS[names(EXPDESIGNOBJ$expDesign ),]),cex.lab=1, las=2 )
		grid()
	}else{
		plot(0,0,main="C.V. NOT CALCULATED")
	}
	
	sampleIntSumBarplot(UNGROUPEDRAWINTDATA
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGNOBJ$expDesign )
			, main="RAW INT SUM PER SAMPLE")
	
	
	par(mar=c(4.1,4.1,1.1,2.1))
	
	plotIntensityDistributions(log10(UNGROUPEDRAWINTDATA)
			,isLegend=F
			, xlab="RAW INTENSITY"
			, ylab="FREQUENCY"
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGNOBJ$expDesign )
			, lwd=1.5)
	
	
	plotIntensityDistributions(log10(UNGROUPEDNORMINTDATA)
			, isLegend=T
			, naReplacedInt=log10(BASELINEINTENSITY)
			, xlab="NORM INTENSITY"
			, ylab="FREQUENCY"
			, colors=getGroupedSampleColors(CONDITIONCOLORS,expDesign=EXPDESIGNOBJ$expDesign )
			, lwd=2)
	
	par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
	
#	plotIntensityDistributions(log10(MEDIANNORMINTDATA)
#			, isLegend=T
#			, naReplacedInt=log10(BASELINEINTENSITY)
#			, xlab="PROTEIN COND. MEDIAN NORM INTENSITY"
#			, ylab="FREQUENCY"
#			, colors=as.character(CONDITIONCOLORS[,1])
#			, lwd=2)
	
}

### INTENSITY DISTRIBUTION PLOTS END

#### HIERARCHICAL CLUSTERING END

### HIERARCHICAL CLUSTERING ratios
if(userOptions$isHClustPlot){
	
	### if more than 500 features perform clustering on subset
	if(exists("ISPEPTIDEANALYSIS") && ISPEPTIDEANALYSIS){
		
		hClustHeatMap3(UNGROUPEDNORMINTDATA
				,expDesign=EXPDESIGNOBJ$expDesign 
				,conditionColors=CONDITIONCOLORS 
				,controlIntensities=MEDIANNORMINTDATA[,CONTROLCONDITION]
				,main=paste("CLUSTERING OF 500 \nMOST INTENSE FEATURES") 
				,selIndices = order(apply(UNGROUPEDNORMINTDATA,1,median),decreasing=T)[1:(min(c(500,nrow(UNGROUPEDNORMINTDATA))))]
		)
		
	}else{
		
		hClustHeatMap3(UNGROUPEDNORMINTDATA
				,expDesign=EXPDESIGNOBJ$expDesign 
				,conditionColors=CONDITIONCOLORS 
				,controlIntensities=MEDIANNORMINTDATA[,CONTROLCONDITION]
		)
		
	}
	
}
### HIERARCHICAL CLUSTERING END

### VOLCANO PLOTS
if(userOptions$isVolcanoPlots){
	plotAllVolcanoes(RATIOSMEDIANNORMINTPERCOND
			, QVALUESPERCOND
			, CVSPERCOND
			, controlCondition=CONTROLCONDITION
			, nonControlConditions=NONCONTROLCONDITIONS
			, qvalueCutOff=userOptions$deFdrCutoff
			, ratioCutOff=userOptions$ratioCutOff
			, ylab="-log10(q-value)"
	)
}

### for experimental purposes eBayes pvals
if(userOptions$eBayes & userOptions$isVolcanoPlots){
	plotAllVolcanoes(RATIOSMEDIANNORMINTPERCOND
			, PVALUESPERCOND
			, CVSPERCOND
			, controlCondition=CONTROLCONDITION
			, nonControlConditions=NONCONTROLCONDITIONS
			, qvalueCutOff=userOptions$deFdrCutoff
			, ratioCutOff=userOptions$ratioCutOff
			,main="moderated t-test P-VALUES"
			, ylab="-log10(p-value)"
	)
	
}
### VOLCANO PLOTS END

### D.E. FDR PLOTS
if(userOptions$isDeFdrPlot){
	
	### two plots per slide	
	par(mfrow=c(1,2))
	
	plotNbValidDeFeaturesPerFDR(QVALUESPERCOND
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=FALSE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalRange=c(0,0.15)
			,qvalCutOff = userOptions$deFdrCutoff
			,conditionColors= CONDITIONCOLORS	
			, main="DOWN-REG. FEATURES"
			,ylab="# FEATURES"
			,xlab="False Discovery Rate (q-value)"
	
	)
	
	plotNbValidDeFeaturesPerFDR(QVALUESPERCOND
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=TRUE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalRange=c(0,0.15)
			,qvalCutOff = userOptions$deFdrCutoff
			,conditionColors= CONDITIONCOLORS
			, main="UP-REG. FEATURES"
    		,ylab="# FEATURES"
			,xlab="False Discovery Rate (q-value)"
			,isLegend=F
	)
	
	
	par(mfrow=c(1,1))
}

### experimental
if(userOptions$eBayes & userOptions$isDeFdrPlot){
	
	### two plots per slide	
	par(mfrow=c(1,2))
	
	plotNbValidDeFeaturesPerFDR(PVALUESPERCOND
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=FALSE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalRange=c(0,0.15)
			,qvalCutOff = userOptions$deFdrCutoff
			,conditionColors= CONDITIONCOLORS	
			, main="DOWN-REG. FEATURES"
			,ylab="# FEATURES"
			,xlab="P-value \n(moderated t-test)"
	
	)
	
	plotNbValidDeFeaturesPerFDR(PVALUESPERCOND
			,RATIOSMEDIANNORMINTPERCOND
			,upRegulated=TRUE
			,logRatioCufOff=log(userOptions$ratioCutOff)
			,qvalRange=c(0,0.15)
			,qvalCutOff = userOptions$deFdrCutoff
			,conditionColors= CONDITIONCOLORS
			, main="UP-REG. FEATURES"
    		,ylab="# FEATURES"
			,xlab="P-value \n(moderated t-test)"
			,isLegend=F
	)
	
	par(mfrow=c(1,1))
	
	### export p-value histograms
	for(cond in names(PVALUESPERCOND)){
		hist(PVALUESPERCOND[,cond], main=paste(cond,"_vs_",CONTROLCONDITION,sep=""), xlab="P-value", breaks=50)
	}
}

### D.E. FDR PLOTS END

cat("CREATED FILE", userOptions$pdfFilePath,"\n")
graphics.off()
