# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### SET COLOR VECTOR
COLORS <- as.character(c("red"
						,"green"
						,"blue"
						,"cyan"
						,"magenta"
						,"darkgreen"
						,"darkblue"
						,"darkorange"
						,"darkgrey"
						,"darkred"
						,"brown"
						,"black"	
						,rev(colors())) ### if that's not enough
		)



### diplay experimental design, high-lighting the control condition 
plotExpDesign <- function(expDesign,sampleNames=sampleNames, controlCondition=controlCondition, condColors=condColors, version="X"){
	
	conditionNames <- names(expDesign)
	nbSamples <- length(sampleNames)
	nbConditions <- length(conditionNames)
	
	xlim <- c(-1,6)
	ylim <- c(-2,NBSAMPLES+2)
	
	plot(0,0,type="n", xlim=xlim, ylim=ylim, main="Experimental Design", axes=FALSE, xlab="", ylab="")
	
	condYPosStep <- (nbSamples+2)/(nbConditions+1)
	sampleNb <- 1
	
	for(condNb in 1:length(conditionNames)){
		
		condName <- conditionNames[condNb]
		condCol = as.character(condColors[condName,])
		
		#cond <- names(EXPDESIGN)[condName]
		### control condition in  and underline
		if(condName == controlCondition){
			text(1,(condNb)*condYPosStep,bquote(underline(bold(.(condName)))), col=condCol)
		}else{
			text(1,(condNb)*condYPosStep,condName, col=condCol)
			
		}
		
		for(i in 1:expDesign[[condName]]){
			text(4,sampleNb,sampleNames[sampleNb], col=condCol)
			sampleNb <- sampleNb + 1
		}
	}
	
	text(-1,-2,paste("SafeQuant v.", version), pos=4)
} 

### idnetification score vs. FDR

plotIdScoreVsFDR <-function(scores,qvals, ...){
	plot(sort(scores),rev(sort(qvals)), ... )
}


### plot identification scores ROC-curve, fdr vs. # valid identifications  
plotROC <- function(qvals,ylab="Valid Id's", fdrMax=1,nbDataPoints=100 ){
	
	fdrRange = c(0,fdrMax)
	fdrCutoffs = seq(0,fdrMax,fdrMax/nbDataPoints)
	
	validIds = c()
	for(fdr in fdrCutoffs){
		
		validIds = c(validIds,sum(qvals < fdr))
		
	}
	
	title <- paste(sum(qvals < 0.01),"FDR 0.01")
	
	plot(fdrCutoffs,sort(validIds), ylab=ylab,xlab="FDR",main=title, xlim=fdrRange,type="l")
	
}

# plot identifications target decoy distribution
plotScoreDistrib <-function(targetScores,decoyScores,nbBins=100,scoreName="score",ylab="counts",title="Score Distribution", xlim=c(0,max(targetScores))){
	
	if(length(targetScores) >0 & length(decoyScores) >0){
		
		minScore = min(c(targetScores,decoyScores))
		maxScore = max(c(targetScores,decoyScores))
		
		binWidth = (maxScore - minScore)/nbBins  
		breaks = seq(minScore - binWidth, maxScore+binWidth,binWidth)
		
		targetHist = hist(targetScores, breaks=breaks, plot=FALSE)
		decoyHist = hist(decoyScores, breaks=breaks, plot=FALSE)
		
		ylim = c(0,max(c(decoyHist$counts,targetHist$counts)))
		
		plot(0,0,type='n', ylim=ylim, xlab=scoreName, ylab=ylab, main=title, xlim=xlim)
		lines(targetHist$mids, targetHist$counts, col=1)
		lines(decoyHist$mids, decoyHist$counts, col=2)
		legend("topright",c("target","decoy"),fill=c(1,2))
	}else{
		print("Not enough target or decoy scores")
	}
}


### plot lower triangle Pearsons R^2. Diagonal text, upper triangle all against all scatter plots with lm abline
pairs.annot<-
		function(data, main="", diagText=c(),...) {
	
	### re-format column names
	names(data) <- gsub("normInt_","", names(data))
	names(data) <- gsub("\\.[1-9]*$","", names(data), perl=T)
	
	count <- 1
	
	panel.lm <-
			function (x, y, col = par("col"), bg = NA, pch = par("pch"),
					cex = 1, col.lm = "red", lwd=par("lwd"), ...)
	{
		#points(x, y, pch = pch, col = col, bg = bg, cex = cex)
		
		points(x, y, pch = 20, col = rgb(0,100,0,50,maxColorValue=255), bg = bg, cex = cex)
		
		
		ok = is.finite(x) & is.finite(y)
		if (any(ok)){
			abline(lm(y~x,subset=ok), col = "red", lwd=1.5, ...)
			abline(coef=c(0,1),lty=2)
		}
	}
	
	panel.sse <-
			function(y, x, digits=2)
	{
		
		usr = par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		
		model = summary(lm(y~x))
		r2= model$r.squared

		txt = round(r2, digits)
		txt = bquote(R^2*"" == .(txt))
		
		text(0.5, 0.5, txt, cex=1.8)
	
	}
	
	panel.txt <- function(x, y, labels, cex, font,digits=2, ...){
		
		txt = colnames(data)[count]
		text(0.5, 0.6, txt, cex=1.8)
		
		if(length(diagText) > 0){
			txt = round(diagText[count], digits)
			txt = bquote(gf == .(txt))
			text(0.5, 0.4, txt, cex=2)
		}
		count <<-count+1
	}
	
	pairs(data,lower.panel=panel.sse,upper.panel=panel.lm, text.panel=panel.txt,main=main, ...)
}


### plot intensity distributions
plotIntensityDistributions <- function(df, isLegend=T, naReplacedInt=-1, colors=colors, ... ){
	
	names= names(df)
	
	maxIntensity <- max(df[df > 0], na.rm=TRUE)[1] +1
	minIntensity <- min(df[df > 0], na.rm=TRUE)[1] -1
	
	histograms <- list()
	
	### get plots y-max (mhc)
	maxHistCounts <- 0
	breaks <- seq(minIntensity,maxIntensity,0.4)
	for(i in 1:length(names)){
		
		data <- df[,i][df[,i] > 0]
		h <- hist(data,plot=FALSE, breaks=breaks)
		histograms[[i]] <- h
		
		mhc <- max(h$counts)[1]
		if(mhc > maxHistCounts){
			maxHistCounts <- mhc 
		}
	}
	
	### plot frame by limits maxIntensity, minIntensity, mhc 
	plot(histograms[[1]]$mids,histograms[[1]]$counts
			,type="n"
			,xlim=c(minIntensity,maxIntensity)
			,ylim=c(0,maxHistCounts),
			...
	)
	
	### draw histogram contours
	for(i in 1:length(names)){
		lines(histograms[[i]]$mids,histograms[[i]]$counts ,col=colors[i], ...)
	}
	
	### draw minIntensity flag
	if(naReplacedInt >= 0){
		
		naFraction <- sum(df[,1] == naReplacedInt) / length(df[,1]) 
		textYPos <- mhc* naFraction * 8
		
		points(naReplacedInt,textYPos, type="h" ,lty=3)
		text(naReplacedInt,textYPos,"Replaced\nNA's", pos=3,cex=0.9)
		
	}
	
	### add legend
	if(isLegend){
		legend("topleft", names ,fill=colors, cex=0.8, pt.cex=0.8 )
	}
}

### barplot of intensity sums per condition
sampleIntSumBarplot <- function(intDf, colors=colors, ...){
	
	sampleNames <- names(intDf)
	sampleNames <- gsub("\\.[0-9]*$","",sampleNames,perl=T)
	
	mp <- barplot(apply(intDf,2,FUN=function(t){return(sum(t,na.rm=TRUE) )}),
			, xaxt = "n"
			, las=2
			, col=colors
			, ...
	)
	
	### 45 degree labels
	axis(1, labels = FALSE)
	text(mp , par("usr")[3] -0.5, srt = 45, adj = 1,
			labels = sampleNames, xpd = TRUE)	
	
}

### get each condition color repeated by number of samples 
getGroupedSampleColors <- function(condColors,expDesign=expDesign){
	
	colors <- c()
	for(cond in names(expDesign)){
		colors <- c(colors,as.character(rep(condColors[cond,],expDesign[[cond]])))
	}
	
	return(colors)
	
}

### color strip for volcano plot
cvColorstrip <- function(colors, maxCV = maxCV  )
{
	BOTTOM <- 1	
	count <- length(colors)
	m <- matrix(1:count, count, 1)
	image(m, col=colors, ylab="", axes=FALSE)
	
	labels <- round(seq(0,maxCV,length.out=3))
	
	### round off to nearest 5 
	labels <- as.character( round(labels / 5) *5 )
	
	at <- seq(0,1,length.out=3)
	
	axis(BOTTOM,tick=TRUE, labels=labels , at=at )
	mtext("C.V. (100%)", BOTTOM, adj=0.5, line=2)
}


### plot log2 ratio vs. abs log10 qvalue, color points by max per condition pair c.v.
volcanoPlot <- function(ratios,qvalues,cvMaxPerPair,cvColorPalette=cvColorPalette,cond=cond ,controlCond=controlCond,qvalueCutOff=qvalueCutOff,ratioCutOff=ratioCutOff  ){
	
	### plot extras params
	ratioCutOffAbsLog2 <- abs(log2(ratioCutOff))
	CVcolors <- cvColorPalette[round(1+ cvMaxPerPair/max(cvMaxPerPair) *31)] 
	
	qvaluesAbsLog <- abs(log10(qvalues)) 
	ratiosLog <- log2(ratios)
	
	par(fig=c(0,1,0.18,1))
	plot(ratiosLog,qvaluesAbsLog,xlab= paste("log2(",cond,"/",controlCond,")", sep="" ), ylab="log10(q-value)", type="n")
	
	points(ratiosLog,qvaluesAbsLog, col=CVcolors, pch=20,)
	validSel <- (qvalues < qvalueCutOff) & (abs(ratiosLog) > ratioCutOffAbsLog2 )
	#points(ratiosLog[validSel],qvaluesAbsLog[validSel], col="lightgrey", pch=1)
	points(ratiosLog[validSel],qvaluesAbsLog[validSel], col=CVcolors, pch=20)
	
	### draw valid squares contianing valid proteins/peptides 
	lines(c(-ratioCutOffAbsLog2,-ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),1000), col="grey")
	lines(c(ratioCutOffAbsLog2,ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),1000), col="grey")
	lines(c(-1000,-ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),abs(log10(qvalueCutOff))), col="grey")
	lines(c(ratioCutOffAbsLog2,1000),c(abs(log10(qvalueCutOff)),abs(log10(qvalueCutOff))), col="grey")
	
	### add heat color bar
	par(fig=c(0,1,0,0.3), new=TRUE)
	cvColorstrip(cvColorPalette, maxCV=max(cvMaxPerPair)*100)
	par(mfrow=c(1,1))
}

### callvolcanoPlot for all vs. control condition pairs
plotAllVolcanoes <- function(ratiosMedianNormIntPerCond,qvaluesPerCond,cvsPerCond,controlCondition=controlCondition,nonControlConditions=nonControlConditions, qvalueCutOff =0.01, ratioCutOff=2){
	
	cvColorPalette <- rev(rich.colors(32))
	
	### iterate over all nonControlConditions
	for(cond in nonControlConditions){
		
		qvalues <- qvaluesPerCond[,cond]
		ratios <- ratiosMedianNormIntPerCond[,cond]
		cvMaxPerPair <-  apply(cvsPerCond[, names(cvsPerCond)  %in% c(cond,controlCondition) ],1,FUN=function(cvs){return(max(cvs, na.rm=T))})
		
		volcanoPlot(ratios,qvalues,cvMaxPerPair,cvColorPalette=cvColorPalette,cond=cond ,controlCond=controlCondition, qvalueCutOff = qvalueCutOff,ratioCutOff=ratioCutOff )
		
	}
}

### hierc clustering heat map normalizing intensities protein wise, sample with highiesintensity assigned intensity 1
hClustHeatMap <- function(data, expDesign=expDesign , conditionColors=conditionColors, ...){
	
	### tidy up sample names
	names(data)  <- gsub("\\.[0-9]*$","",names(data), perl=TRUE )
	names(data)  <- gsub("normInt\\_","",names(data), perl=TRUE )
	
	### NORMALIZE PER PROTEIN
	proteinMaxIntensity = apply(data, 1 ,FUN=max)
	normPerProteinData = data / proteinMaxIntensity
	
	### BI-CLUSTERING HEAT MAP
	protein.cor = cor(t(normPerProteinData), use="pairwise.complete.obs", method="pearson")
	protein.cor.dist = as.dist(1-protein.cor)
	proteins.tree = hclust(protein.cor.dist, method='average')
	
	msrun.cor.spearman = cor(normPerProteinData, use="pairwise.complete.obs", method="spearman")
	msrun.cor.spearman.dist = as.dist(1-msrun.cor.spearman)
	msrun.tree = hclust(msrun.cor.spearman.dist,method='average')
	
	### sample colors
	samplecolors = getGroupedSampleColors(conditionColors,expDesign=expDesign)
	
	labRow <- rownames(data)
	
	### do not display feature names if too many
	if(nrow(data) > 50){
		labRow <- rep("",(nrow(data)))
	}
	
	heatmap.2(as.matrix(normPerProteinData)
			, col=colorRampPalette(c(colors()[142],"black",colors()[128]))
			, scale="none"
			, key=TRUE
			, symkey=FALSE
			#, density.info="none"
			, trace="none"
			, cexRow=0.5
			, cexCol = 0.7
			,ColSideColors=samplecolors
			,labRow = labRow
			,Rowv=as.dendrogram(proteins.tree)
			,Colv=as.dendrogram(msrun.tree)
			, ...
	)
	
	legend("topright", names(expDesign), fill=as.character(conditionColors[,1]), cex=0.6)
}

### plot qvalue vs. tot. valid features (ROC-curver), applying ratio cutoff. plotting all conditions
plotNbValidDeFeaturesPerFDR <- function(qvaluesPerCond,ratiosPerCond,upRegulated=T,logRatioCufOff=log(1.5),qvalCutOffs=seq(0,0.3, length.out=10),conditionColors= conditionColors, isLegend=T, ... ){
	
	conditions <- names(qvaluesPerCond)
	
	### create data farme of roc curve per condition
	### store curves
	nbPassingCutOffsPerCond <- data.frame(row.names=qvalCutOffs)	
	
	### iterate over all conditions
	for(cond in conditions){
		
		qvals <- qvaluesPerCond[cond]	
		ratios <- ratiosPerCond[cond]
		ratiosLog <- log(ratios)
		
		#iterate over all cutoffs
		nbPassingCutOffs <- c()
		for(qCutOff in qvalCutOffs){
			
			if(upRegulated){
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( qvals < qCutOff & ratiosLog >  logRatioCufOff ) )
			}else{
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( qvals < qCutOff & ratiosLog <  -logRatioCufOff ) )
			}
		}
		
		nbPassingCutOffsPerCond <- data.frame(nbPassingCutOffsPerCond, nbPassingCutOffs)
		
	}
	
	names(nbPassingCutOffsPerCond) <- conditions
	
	# plot roc curves	
	plot(0,0, ylim=c(0,max(nbPassingCutOffsPerCond)), xlim= c(min(qvalCutOffs) , max(qvalCutOffs)), type="n", xlab="False Discovery Rate (q-value)", ...)
	for(cond in conditions){
		lines(qvalCutOffs,nbPassingCutOffsPerCond[,cond], col= as.character(conditionColors[cond ,]), lwd=1.6)
	}
	if(isLegend){
		legend("topleft", conditions, fill=as.character(conditionColors[conditions,]), cex=0.6)
	}
}


### plot mass diff v.s score highlighting decoy and displaying user specified user specified precursor mass filter
plotMassDiffVsScore <- function(massDiff, scores, decoyCond, precursorMassFilter=precursorMassFilter ){
	
	plot(massDiff, scores 
			, col = decoyCond+1
			, pch=20, ylab="score", xlab="mass diff. (ppm)")
	
	lines(c(precursorMassFilter,precursorMassFilter),c(min(scores),max(scores)),col="blue",lwd=2)
	legend("topright",c("target","decoy"), pch=20, col= c(1,2))
	
}