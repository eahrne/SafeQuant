# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### SET COLOR VECTOR
COLORS <- as.character(c(
				"red"
				,"darkgreen"
				,"blue"
				,"darkmagenta"
				,"darkorange"
				,"darkblue"
				,"black"
				,"brown"
				,"darkgrey"
				,"red"
				,"brown"
				,rev(colors())) ### if that's not enough
)

### some quality control plots
.qcPlots <- function(eset,selection=1:6,nbFeatures=500, ...){
	
	if(nrow(eset) < nbFeatures) nbFeatures <- nrow(eset) 
	
	sel <- sample(nrow(eset),nbFeatures,replace=F)
	
	if(1 %in% selection) plotExpDesign(eset)
	if(2 %in% selection)barplotMSSignal(eset, ...)
	if(3 %in% selection)plotMSSignalDistributions(log2(exprs(eset)), col=as.character(.getConditionColors(eset)[pData(eset)$condition,]), lwd=1.5, ...)
	if(4 %in% selection)pairsAnnot(log2(exprs(eset)[sel,order(pData(eset)$condition)]), ...)
	if(5 %in% selection)hClustHeatMap(eset[sel,], ...)
	### retention time normalization plot
	if((6 %in% selection) && ("rt" %in% sqaMethod) ){
		plotRTNormSummary(getRTNormFactors(eset, minFeaturesPerBin=100))
	}
	
	
} 


.getConditionColors <- function(eset){
	return(data.frame(colors=as.character(COLORS[1:length(levels(pData(eset)$condition))]), row.names=levels(pData(eset)$condition)))	
}

### color strip for volcano plot
.cvColorstrip <- function(colors, maxCV = maxCV  )
{
	bottom <- 1	
	count <- length(colors)
	m <- matrix(1:count, count, 1)
	image(m, col=colors, ylab="", axes=FALSE)
	
	labels <- round(seq(0,maxCV,length.out=3))
	
	### round off to nearest 5 
	labels <- as.character( round(labels / 5) *5 )
	
	at <- seq(0,1,length.out=3)
	
	axis(bottom,tick=TRUE, labels=labels , at=at )
	mtext("C.V. (%)", bottom, adj=0.5, line=2)
}

### called from plotVolcano. Creates volcano plot form data.frame input
.plotVolcano <- function(d
		, ratioCutOffAbsLog2=0
		, absLog10pValueCutOff=2
		, xlim = range(d[,1],na.rm=T)
		, ylim = range(abs(        log10(d[,2])[is.finite(  log10(d[,2])  )]    )      ,na.rm=T) # pValue or qValue
		, cvMax = max(d[,3][is.finite(d[,3])],na.rm=T)	
		#, higlightSel = rep(F,nrow(d))
		, controlCondition = "control"
		, caseCondition = "case"
		,...

){
	
	cvColorPalette <- rev(rich.colors(32))
	CVcolors <- cvColorPalette[round(1+ d[,3]/cvMax *31)] 
	
	par(fig=c(0,1,0.18,1))
	plot(0,0
			, xlab= paste("log2(",caseCondition,"/",controlCondition,")", sep="" )
			, ylab= paste("-log10(",names(d)[2],")",sep="")
			, type="n", ylim=ylim,xlim=xlim	
			,...
	)
	
	grid()
	
	### d[,2] qValues or pValues
	points(d[,1], abs(log10(d[,2])), col=CVcolors, pch=20)
	
	### draw valid squares contianing valid proteins/peptides 
	lines(c(-ratioCutOffAbsLog2,-ratioCutOffAbsLog2),c(absLog10pValueCutOff,1000), col="grey")
	lines(c(ratioCutOffAbsLog2,ratioCutOffAbsLog2),c(absLog10pValueCutOff,1000), col="grey")
	lines(c(-1000,-ratioCutOffAbsLog2),c(absLog10pValueCutOff,absLog10pValueCutOff), col="grey")
	lines(c(ratioCutOffAbsLog2,1000),c(absLog10pValueCutOff,absLog10pValueCutOff), col="grey")
	
	# @TODO remove this feature
	#points(d$ratio[higlightSel],abs(log10(d[,2]))[higlightSel],col="darkgrey")
	
	### add heat color bar
	par(fig=c(0,1,0,0.3), new=TRUE)
	.cvColorstrip(cvColorPalette, maxCV=cvMax*100)
	par(mfrow=c(1,1))
	
}

#import gplots
#' Plots volcano, data points colored by max cv of the 2 compared conditions
#' @param obj safeQuantAnalysis object or data.frame
#' @param ratioCutOffAbsLog2 ratio abline 
#' @param absLog10pValueCutOff pValue abline
#' @param  adjusted TRUE/FALSE plot qValues or pValues on y-axis
#' @export
#' @import Biobase gplots
#' @note  No note
#' @details data.frame input object should contain 3 columns (ratio,qValue,cv)
#' @references NA
#' @examples print("No examples")
plotVolcano <- function(obj
		,ratioCutOffAbsLog2=0
		,absLog10pValueCutOff=2
		,adjusted = T
		
		,...
){
	### plot volcanos for all case control comparisons
	if(class(obj) == "safeQuantAnalysis"){
		
		# ensure the same range on all volcano plots
		xlim <- range(obj$ratio, na.rm=T)
		if(adjusted){
			ylim <- range(abs(log10(obj$qValue)),na.rm=T)
		}else{
			ylim <- range(abs(log10(obj$pValue)),na.rm=T)
		}
		# ensure same scale on color legend
		cvMax <- max(obj$cv,na.rm=T)
		controlCondition <- .getControlCondition(obj$eset)
		
		for(caseCondition in colnames(obj$pValue)){
			
			### create data.frame listing all data points (ratio,pvalue,cv)
			if(adjusted){
				
				d <- data.frame( ratio=as.vector(unlist(obj$ratio[caseCondition]))
						,qValue=as.vector(unlist(obj$qValue[caseCondition]))
						,cv=apply(obj$cv[c(caseCondition,controlCondition)],1,max,na.rm=T)
				)
			}else{
				d <- data.frame( ratio=as.vector(unlist(obj$ratio[caseCondition]))
						,pValue=as.vector(unlist(obj$pValue[caseCondition]))
						,cv=apply(obj$cv[c(caseCondition,controlCondition)],1,max,na.rm=T)
				)
			}
			
			.plotVolcano(d
					, ratioCutOffAbsLog2=ratioCutOffAbsLog2
					, absLog10pValueCutOff=absLog10pValueCutOff
					, xlim = xlim
					, ylim = ylim # pValue or qValue
					, cvMax = cvMax	
					#, higlightSel = higlightSel
					, controlCondition = controlCondition
					, caseCondition = caseCondition
					
					,...)
		} 
	}else{
		.plotVolcano(obj
				, ratioCutOffAbsLog2=ratioCutOffAbsLog2
				, absLog10pValueCutOff=absLog10pValueCutOff
				#, higlightSel = higlightSel
				
				,...)
	}
	
}

#' Display experimental design, high-lighting the control condition 
#' @param eset ExpressionSet
#' @export
#' @import Biobase 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotExpDesign <- function(eset, condColors=.getConditionColors(eset),  version="X"){
	
	conditionNames <- as.character(unique(pData(eset)$condition))
	nbConditions <- length(conditionNames)
	controlCondition = .getControlCondition(eset)
	
	sampleNames <- rownames(pData(eset))
	nbSamples <- nrow(pData(eset))
	
	xlim <- c(-1,6)
	ylim <- c(-2,nbSamples+2)
	
	plot(0,0,type="n", xlim=xlim, ylim=ylim, main="Experimental Design", axes=FALSE, xlab="", ylab="")
	
	condYPosStep <- (nbSamples+2)/(nbConditions+1)
	sampleNb <- 1
	
	for(condNb in 1:length(conditionNames)){
		
		condName <- conditionNames[condNb]
		condCol = as.character(condColors[condName,])
		
		### control condition in  and underline
		if(condName == controlCondition){
			text(1,(condNb)*condYPosStep,bquote(underline(bold(.(condName)))), col=condCol)
		}else{
			text(1,(condNb)*condYPosStep,condName, col=condCol)
			
		}
		
		for(i in sampleNames[as.character(pData(eset)$condition) == condName]){
			text(4,sampleNb,paste(paste(sampleNb,":",sep=""), sampleNames[sampleNb]), col=condCol)
			sampleNb <- sampleNb + 1
		}
	}
	
	text(-1,-2,paste("SafeQuant v.", version), pos=4)
	
}

#' Plot lower triangle Pearsons R^2. Diagonal text, upper triangle all against all scatter plots with lm abline
#' @param data data.frame
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
pairsAnnot<-
		function(data, diagText=c(),...) {
	
	if(ncol(data) < 2 ){
		cat("ERROR: pairsAnnot, min 2 data columns required \n")
		return(-1)
	}
	
	### re-format column names
	names(data) <- gsub("normInt_","", names(data))
	### hack
#	if( regexpr("Condition", names(data)[1]) == -1 ){
#		names(data) <- gsub("\\.[1-9]{1,2}$","", names(data), perl=T)
#	}
#	
	count <- 1
	
	panel.lm <-
			function (x, y, col = par("col"), bg = NA, pch = par("pch"),
					#	cex = 1, col.lm = "red", lwd=par("lwd"), ...)
					col.lm = "red", lwd=par("lwd"))
	{
		points(x, y, pch = 20, col = rgb(0,100,0,50,maxColorValue=255), bg = bg)
		
		ok = is.finite(x) & is.finite(y)
		if (any(ok)){
			abline(lm(y~x,subset=ok), col = "red", lwd=1.5)
			abline(coef=c(0,1),lty=2)
		}
	}
	
	panel.sse <-
			function(y, x, digits=2,...)
	{
		
		usr = par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		
		### discard non-fininte values
		ok = is.finite(x) & is.finite(y)
		x <- x[ok]
		y <- y[ok]
		
		model = summary(lm(y~x))
		r2= model$r.squared
		
		txt = round(r2, digits)
		txt = bquote(R^2*"" == .(txt))
		
		text(0.5, 0.5, txt)
		
	}
	
	#panel.txt <- function(x, y, labels, cex, font,digits=2, ...){
	panel.txt <- function(x, y, labels, font,digits=2,...){
		
		txt = colnames(data)[count]
		#text(0.5, 0.6, txt, cex=1.5)
		text(0.5, 0.6, txt)
		
		if(length(diagText) > 0){
			txt = round(diagText[count], digits)
			txt = bquote(gf == .(txt))
			text(0.5, 0.4, txt)
			#text(0.5, 0.4, txt, cex=1.5)
		}
		count <<-count+1
	}
	
	pairs(data,lower.panel=panel.sse,upper.panel=panel.lm, text.panel=panel.txt,...)
}

#' Plot ms.signal distributions
#' @param matrix matrix of ms-signals
#' @param color color
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotMSSignalDistributions <- function(d, col=1:100, cex.axis=1, cex.lab=1,ylab="binned count", xlab="MS-signal",... ){
	
	
	
	### create a sms-signal histogram per column,
	breaks <- seq(min(d,na.rm=T),max(d,na.rm=T),length=40)
	mids <-  hist(d[,1],breaks=breaks,plot=F)$mids
	countPerCol <- data.frame(row.names=mids)
	ncol(countPerCol)
	
	for(c in colnames(d)){
		count <- hist(d[,c],breaks=breaks,plot=F)$count
		countPerCol <- cbind(countPerCol,count)
	}
	names(countPerCol) <-  colnames(d)
	
	### plot histogram trend lines
	plot(mids,mids,ylim=range(countPerCol)
			, type="n"
			, cex.axis=cex.axis
			, cex.lab=cex.lab
			,xlab=xlab
			,ylab=ylab
			,...)
	for(i in 1:ncol(countPerCol)){
		lines(mids,countPerCol[,i],col=col[i],...)
	}
	legend("topleft", names(countPerCol),lty=1 , col=col)
	
}

### 
#' Barplot of ms-signal sums per column
#' @param matrix matrix of ms-signals
#' @param color color
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
barplotMSSignal <- function(eset, cex.lab=1,...){
	
	### handle 0 and negavtive values
	exprs(eset)[exprs(eset) <= 0] <- NA
	
	#d <- log10(exprs(eset))
	colors<-as.vector(unlist(.getConditionColors(eset))[pData(eset)$condition])
	
	mp <- barplot( apply(exprs(eset),2,FUN=function(t){return(sum(t,na.rm=TRUE) )}),
			, xaxt = "n"
			, las=2
			, col=colors
			,names=names(exprs(eset))
			, ...
	)
	
	### 45 degree labels
	axis(1, labels = FALSE,tick=F)
	#text(mp , par("usr")[3] -0.5, srt = 45, adj = 1.2,
	text(mp , par("usr")[3], srt=45, adj = 1,
			labels = colnames(exprs(eset)), xpd = TRUE, ,cex=cex.lab)	
	
}

### 
#' Hierc. clustering heat map, cluster by mms signal, display log2 ratios to control median
#' @param eset ExpressionSet
#' @param conditionColors data.frame of colors per condition
#' @param selIndices indices of selected 
#' @return heatmap.2 obj
#' @export
#' @import ExpressionSet gplots
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
hClustHeatMap <- function(eset
		,conditionColors =.getConditionColors(eset)
		,selIndices = 1:nrow(eset)
		,breaks=seq(-2,2,length=20)
		,...
){
	
	d <- log2(exprs(eset)[selIndices,])
	### BI-CLUSTERING HEAT MAP
	
	feature.cor = cor(t(d), use="pairwise.complete.obs", method="pearson")
	feature.cor.dist = as.dist(1-feature.cor)
	feature.cor.dist[is.na(feature.cor.dist)] <- 0
	feature.tree = hclust(feature.cor.dist, method='median')
	
	msrun.cor.pearson = cor(d, use="pairwise.complete.obs", method="pearson")
	msrun.cor.pearson.dist = as.dist(1-msrun.cor.pearson)
	### to avoid error when replicates of the same condition are identical, DOES THIS EVER HAPPEN?
	msrun.cor.pearson.dist[is.na(msrun.cor.pearson.dist)] <- 0
	msrun.tree = hclust(msrun.cor.pearson.dist,method='median')
	
	### sample colors
	samplecolors =  as.vector(unlist(conditionColors[pData(eset)$condition,]))
	
	### log2 ratios to median of control condition
	log2RatioPerMsRun <- d - log2(getSignalPerCondition(eset,method="median")[,.getControlCondition(eset)])
	
	labRow <- rownames(log2RatioPerMsRun)
	
	### do not display feature names if too many
	if(nrow(log2RatioPerMsRun) > 50){
		labRow <- rep("",(nrow(log2RatioPerMsRun)))
	}
	
	hm <- heatmap.2(as.matrix(log2RatioPerMsRun)
			, col=colorRampPalette(c(colors()[142],"black",colors()[128]))
			, scale="none"
			, key=TRUE
			, symkey=FALSE
			, trace="none"
			, cexRow=0.5
			, cexCol = 0.7
			,ColSideColors=samplecolors
			,labRow = labRow
			,Rowv=as.dendrogram(feature.tree)
			,Colv=as.dendrogram(msrun.tree)
			,dendrogram="column"
			,density.info="density"
			#,KeyValueName="Prob. Response"
			,breaks=breaks
			, ...
	)
	
	legend("left",levels(pData(eset)$condition), fill=as.character(conditionColors[,1]), cex=0.7, box.col=0)
	
	return(hm)
	
}



### 
#' Plot Total Number of diffrentially Abundant Features (applying ratio cutoff) vs. qValue/pValue for all conditions
#' @param sqa SafeQuantAnalysis Object
#' @param upRegulated TRUE/FALSE select for upregulated features 
#' @param log2RatioCufOff log2 ratio cut-off
#' @param pvalCutOff pValue/qValue cut-off 
#' @param pvalRange pValue/qValue range
#' @param isLegend TRUE/FALSE display legend
#' @param isAdjusted TRUE/FALSE qValues/pValue on x-axis
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotNbValidDeFeaturesPerFDR <- function(sqa,upRegulated=T,log2RatioCufOff=log2(1.5),pvalRange=c(0,0.3)
	,pvalCutOff=1, isLegend=T,isAdjusted=T,ylab="Nb. Features", ... ){
	
	if(isAdjusted){
		pvaluesPerCond <- sqa$qValue
		xlab= "qValue"
	}else{
		pvaluesPerCond <- sqa$pValue
		xlab= "pValue"
	}

	ratiosPerCond <- sqa$ratio
	conditionColors <- .getConditionColors(sqa$eset)

	pvalCutOffs <- seq(pvalRange[1],pvalRange[2], length.out=10)
	conditions <- names(pvaluesPerCond)
	
	### create data farme of roc curve per condition
	### store curves
	nbPassingCutOffsPerCond <- data.frame(row.names=pvalCutOffs)	
	
	### iterate over all conditions
	for(cond in conditions){
		
		pvals <- pvaluesPerCond[cond]	
		ratios <- ratiosPerCond[cond]
			
		#iterate over all cutoffs
		nbPassingCutOffs <- c()
		for(qCutOff in pvalCutOffs){
			
			if(upRegulated){
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( (pvals < qCutOff) & (ratios >  log2RatioCufOff) ,na.rm=T) )
			}else{
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( (pvals < qCutOff) & (ratios <  -log2RatioCufOff) ,na.rm=T ) )
			}
		}
		
		nbPassingCutOffsPerCond <- data.frame(nbPassingCutOffsPerCond, nbPassingCutOffs)
		
	}
	
	names(nbPassingCutOffsPerCond) <- conditions
	
	# plot roc curves	
	plot(0,0, ylim=c(0,max(nbPassingCutOffsPerCond)), xlim= c(min(pvalCutOffs) , max(pvalCutOffs)), type="n",ylab=ylab,xlab=xlab,  ...)
	grid()
	for(cond in conditions){
		lines(pvalCutOffs,nbPassingCutOffsPerCond[,cond], col= as.character(conditionColors[cond ,]), lwd=1.6)
	}
	if(isLegend){
		legend("topleft", conditions, fill=as.character(conditionColors[conditions,]), cex=0.6)
	}
	
	abline(v=pvalCutOff,col="grey",lwd=1.5)
}




#' Plot identifications target decoy distribution
#' @param targetScores
#' @param decoyScores
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotScoreDistrib <-function(targetScores,decoyScores,xlab="Identification Score",ylab="Counts", ...){
	
	if(length(targetScores) >0 & length(decoyScores) >0){
		
		breaks <- hist(c(targetScores,decoyScores),plot=FALSE,breaks=100)$breaks
		
		targetHist = hist(targetScores, breaks=breaks, plot=FALSE)
		decoyHist = hist(decoyScores, breaks=breaks, plot=FALSE)
		
		ylim = c(0,max(c(decoyHist$counts,targetHist$counts)))
		
		plot(0,0,type='n', ylim=ylim, xlab=xlab, ylab=ylab, xlim=range(targetScores,decoyScores),...)
		grid()
		points(targetHist$mids[targetHist$counts > 0], targetHist$counts[targetHist$counts > 0], col=1, type="h", lwd=4)
		points(decoyHist$mids[decoyHist$counts > 0], decoyHist$counts[decoyHist$counts > 0], col=2, type="h", lwd=5)
		
	}else{
		
		if(length(targetScores) >0){
			targetHist = hist(targetScores, breaks=100, plot=FALSE)
			plot(targetHist$mids, targetHist$counts, xlab=xlab, ylab=ylab,type="n",... )
			points(targetHist$mids[targetHist$counts > 0], targetHist$counts[targetHist$counts > 0], col=1, type="h", lwd=4)
		}
		
		cat("plotScoreDistrib: Not enough target or decoy scores\n")
	}
	
	legend("topright",c("target","decoy"),fill=c(1,2))
}



#' Plot FDR vs. identification score
#' @param idScore vector of identification scores 
#' @param qvals vector of q-valres
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotIdScoreVsFDR <-function(idScore,qvals, ylab="False Discovery Rate", xlab="Identification Score",...){
	plot(sort(idScore),rev(sort(qvals)),type="l",ylab=ylab,xlab=xlab, ... )
	grid()
}



### plot identification scores ROC-curve, fdr vs. # valid identifications 
#' Plot Number of Identifications vs. FDR 
#' @param qvals vector of q-values
#' @param breaks see breaks for hist function
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotROC <- function(qvals
		,xlab="False Discovery Rate"
		,ylab="# Valid Identifications"
		,xlim=c(0,1)
		,breaks=100
		,col="blue"
		,lwd=1.5	
		,... ){
	
	if(length(breaks) == 1 ){
		breaks = seq(xlim[1],xlim[2],length=breaks)
	}else{  ### breaks override xlim
		xlim <- c(min(breaks),max(breaks))
	}

	validIds = c()
	for(fdr in breaks){
		validIds = c(validIds,sum(qvals < fdr))
	}
	
	plot(breaks,sort(validIds), ylab=ylab,xlab=xlab, xlim=xlim,type="l",col=col,lwd=lwd,...)
	grid()
}



### Plot Precursor Mass Error Distribution 
#' Plot Precursor Mass Error Distribution 
#' @param pMassError vector of precursor mass errors 
#' @param isDecoy vector TRUE/FALSE 
#' @param pMassTolWindow Precursor Mass Error Tolerance Window
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotPrecMassErrorDistrib <- function(pMassError,isDecoy,pMassTolWindow=c(-10,10)){
	
	par(mar=c(5.1,4.1,1.1,4.1))
	
	xRange 	<- range(pMassError,pMassTolWindow)
	breaks <- seq(xRange[1],xRange[2],length=50)
	
	### decoy hist
	decoyMdHist 	<- hist(pMassError[isDecoy],breaks=breaks,plot=F)
	
	### target hist
	targetMdHist 	<- hist(pMassError[!isDecoy],breaks=breaks,plot=F)
	
	yRange <- range(decoyMdHist$counts,targetMdHist$counts)
	
	### plot histogram contours
	plot(decoyMdHist$mids,decoyMdHist$counts,type="l", col="red",ylim=yRange,xlab="mass diff. (ppm)",ylab="PSM Frequnecy",lwd=1.5)
	grid()
	lines(decoyMdHist$mids,targetMdHist$counts,lwd=1.5)
	abline(v=pMassTolWindow[1],col="lightgrey",lwd=1)
	abline(v=pMassTolWindow[2],col="lightgrey",lwd=1)
	
	ratio <- (decoyMdHist$counts+1)/(targetMdHist$counts+1)
	ratio[ratio > 1] <- 1
	ratio <- ratio*yRange[2]
	
	lines(decoyMdHist$mids,ratio,lwd=1.5,col="blue",lty=2)
	
	### add ratio ratio line
	axis(4,at=seq(0,yRange[2],length=5), labels=seq(0,1,length=5),col="blue", col.ticks="blue")
	mtext("decoy-target PSM count ratio ",side=4, col="blue",line=2.5)
	
	legend("left",c("target","decoy"),  fill= c("black","red"),bg="white")
	
	par(mar=c(5.1,4.1,4.1,2.1))
	
}

### 
#' Plot precursorMass error v.s score highlighting decoy and displaying user specified user specified precursor mass filter
#' @param pMassError vector of precursor mass errors 
#' @param idScore vector of identification scores
#' @param isDecoy vector TRUE/FALSE 
#' @param pMassTolWindow Precursor Mass Error Tolerance Window
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotPrecMassErrorVsScore <- function(pMassError, idScore, isDecoy, pMassTolWindow=c(-10,10) ){
	
	withinTol <- (pMassError >= pMassTolWindow[1]) & (pMassError <= pMassTolWindow[2])
	
	plot(pMassError, idScore 
			, type="n"		
			#, col = isDecoy+1
			,ylab="score", xlab="mass diff. (ppm)", xlim=range(c(pMassTolWindow,pMassError)))
	
	grid()
	
	points(pMassError[withinTol], idScore[withinTol], col="black", pch=20 )
	points(pMassError[!withinTol], idScore[!withinTol], col="grey", pch=20 )
	points(pMassError[isDecoy], idScore[isDecoy], col="red", pch=20 )
	
	abline(v=pMassTolWindow[1],col="lightgrey",lwd=1)
	abline(v=pMassTolWindow[2],col="lightgrey",lwd=1)
	
	legend("topright",c("target-kept","decoy","target-discarded"), pch=20, col= c("black","red","grey"))
	
}

#' Scatter plot with density coloring
#' @param x number vector
#' @param y number vector
#' @export
#' @import 
#' @note  No note
#' @references NA
#' @examples print("No examples")
plotXYDensity <- function(x,y,isFitLm=T,legendPos="bottomright",disp=c("abline","R"),  ...){
	
	df <- data.frame(x,y)
	
	## Use densCols() output to get density at each point
	densityCol <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
	df$densityCol <- col2rgb(densityCol)[1,] + 1L
	
	## Map densities to colors
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
					"#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$densityCol]
	
	## Plot it, reordering rows so that densest points are plotted on top
	plot(y~x, data=df[order(df$densityCol),], pch=20, col=col, cex=2, ...)
	
	### disp linerar model
	if(isFitLm){
		ok <- is.finite(x) & is.finite(y)
		fit <- lm(y[ok] ~ x[ok])
		
		if("abline" %in% disp){
			abline(coef=c(0,1),lty=2)
			abline(fit)
		}
	
		if("R" %in% disp){
			legend(legendPos
					,legend=c(as.expression(bquote(R^2*"" == .(round(summary(fit)$r.squared,2)))))
					,text.col=c(1,2), box.col=0, cex=2
			)
		}
		return(fit)
	}
	
	return(NA)
	
}

#' Plot calibration Curve
#' @param fit simple log-linear model
#' @param dispElements c("formula","lowess","stats")
#' @export
#' @import 
#' @note  No note
#' @references NA
#' @examples print("No examples")
plotCalibrationCurve <- function(fit
		,dispElements = c("formula","lowess","stats")
		,xlab="Conc. (CPC) "
		,ylab="Pred. Conc."
		,predictorName = paste("log10(",names(coef(fit))[2],")",sep="")
		,text=F,...){
	x <- predict(fit) + fit$residuals 						
	y <- predict(fit)		
	
	
	### some extra margin for axis labels
	par(mar=c(5.5,5.5,4.1,2.1))
	plot(10^x, 10^y,log="xy",xlab="",ylab="",... )
	
	if(text){
		text(10^x, 10^y,rownames(fit$qr$qr))
	}
	
	
	abline(coef=c(0,1),lty=2)
	### add axis labels
	mtext(side=1,xlab,las=1, line=4, ...)
	mtext(side=2,ylab,las=3, line=4, ...)
	
	if( "formula" %in% dispElements){
		legend("bottomright"
				,paste("log10(",ylab,")"," = ", signif(coef(fit)[1],2)," + ",signif(coef(fit)[2],2)," * ",predictorName ,sep="")		
				,box.lwd=0
				,box.col="white"
				,...
		)
	}
	
	if( "lowess" %in% dispElements){
		lws <- lowess(y ~ x)
		lines(10^lws$x, 10^lws$y,col="red",...)
	}
	
	if( "stats" %in% dispElements){
		
		df <- data.frame(cpc =  x,signal = y)
		medianFoldError <- median(abs(getLoocvFoldError(df)[,1]),na.rm=T)
	
		legend("topleft"
				,legend=c(as.expression(bquote(R^2*"" == .(round(summary(fit)$r.squared,2))))
							,paste("Median Fold Error = ",round(medianFoldError,2))
						)
				#,text.col=c(1,2)
				,box.lwd=0
				,box.col="white"
				,...
		)
	}
	### reset margins
	par(mar=c(5.1,4.1,4.1,2.1))
}


#' Plot all retention time normalization profiles
#' @param rtNormFactors data.frame of normalization factor per r.t bin and sample, obtained by getRTNormFactors
#' @param condNames  vector of condition names
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
plotRTNormSummary <- function(rtNormFactors,ylim=range(rtNormFactors), col = 1:ncol(rtNormFactors),condNames=paste("Cond",1:length(unique(col))),...){
	
	plot(rownames(rtNormFactors),rtNormFactors[,1], type="n", ylab="log2(Ratio)",xlab="Retention Time (min)",ylim=ylim, ...)
	for(i in 1:ncol(rtNormFactors)){
		lines(rownames(rtNormFactors),rtNormFactors[,i], col=as.character(col[i]), ...   )
	}
	abline(h=0, lty=2)
	legend("bottom",condNames,lty=1, col=unique(col), lwd=1.5)
}

#' Plot all retention time profile overalying ratios
#' @param rtNormFactors data.frame of normalization factor per r.t bin and sample, obtained by getRTNormFactors
#' @param eset  ExprsssionSet
#' @param samples specify samples (sample numbers) to be plotted
#' @export
#' @import 
#' @note  No note
#' @details No details
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
plotRTNorm <- function(rtNormFactors,eset,samples=1:ncol(rtNormFactors), ...){
	
	### select for anchor proteins
	eset <- eset[fData(eset)$isNormAnchor,]
	
	# get all ratios to sample 1
	# @TODO How to select reference run?
	ratios <- log2(exprs(eset)) - log2(exprs(eset)[,1])
	#ratios <- log2(exprs(eset)) - log2(apply(exprs(eset),1,median,na.rm=T))
	
	for(samplesNb in samples){
		plot(fData(eset)$retentionTime,	ratios[,samplesNb]
				, col="lightgrey"
				, xlab="Retention Time (min)"
				, ylab="log2(Ratio)"
				, main=names(rtNormFactors)[samplesNb], ...)
		lines(as.numeric(rownames(rtNormFactors)),as.vector(unlist(rtNormFactors[,samplesNb])), ... )
		abline(h=0,lty=2,...)
	}
}


## Data in a data.frame

.errorBar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
