# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


### FUNCTIONS

getNormalizationFactors <- function(raw){
	
	#refRunIntSum = sum(raw[,1])
	#rawDataSums = lapply(raw, FUN=sum)
	rawDataSums = apply(raw,2, FUN=sum)
	gainFactors = as.numeric(rawDataSums[1]) / as.numeric(rawDataSums)
	
	return(gainFactors)
	
}

normalizeIntensities <- function(rawData, gainFactors){
	
	normDataColNames = gsub("raw","",names(rawData))
	
	return(matrix(t(apply(rawData, 1, function(.rawData)mapply(gainFactors, .rawData, FUN="*"))),ncol=ncol(rawData),dimnames=list(row.names(rawData),normDataColNames)))
	
}

#normalizeIntensities <- function(rawData){
#	
#	normDataColNames = gsub("raw","",names(rawData))
#	
#	refRunIntSum = sum(rawData[,1])
#	rawDataSums = lapply(rawData, FUN=sum)
#	gainFactors = as.numeric(rawDataSums[1]) / as.numeric(rawDataSums)
#	# multiply each raw data column by corresponding gain factor 
#	normData = matrix(t(apply(rawData, 1, function(.rawData)mapply(gainFactors, .rawData, FUN="*"))),ncol=ncol(rawData),dimnames=list(row.names(rawData),normDataColNames))
#	return(list("normData"=normData, "gainFactors"=gainFactors))
#}

calculateMedianIntensitiesPerCond <- function(e, expDesign){
	
	conditionMedianIntensities <- data.frame(row.names = rownames(e)) 
	minIntensity = min(e)
	
	skip = 1
	for(cond in names(expDesign)){
		
		replicates = expDesign[[cond]]
		#tag <- paste(cond,"MedianIntensity",sep="")
		tag = cond
		xCondFirstReplicate = skip
		xCondLastReplicate = xCondFirstReplicate+replicates-1
		xCondSpectralCount = data[,xCondFirstReplicate:xCondLastReplicate]
		medianIntensity = list()
		#medianIntensity[[tag]] = apply(e[,xCondFirstReplicate:xCondLastReplicate], 1 ,FUN=median)
		
		if(xCondFirstReplicate == xCondLastReplicate){
			
			medianIntensity[[tag]] = e[,xCondFirstReplicate]
			
		}else{
			medianIntensity[[tag]] = apply(e[,xCondFirstReplicate:xCondLastReplicate], 1 ,FUN=function(replicateIntensities){
				
				medianInt = median(replicateIntensities)		
				nbMeasuredIntensities = length(replicateIntensities[replicateIntensities>minIntensity])
				nbReplicates = length(replicateIntensities)
					
				if((nbMeasuredIntensities < nbReplicates) &&  (nbMeasuredIntensities > 0) ){
					medianInt = median(replicateIntensities[replicateIntensities>minIntensity])
	
				}
				return(medianInt)
			})
		}
		conditionMedianIntensities = data.frame(conditionMedianIntensities,medianIntensity)
		skip = skip+replicates
		
	}
	return(conditionMedianIntensities)
}

calculateMeanIntensitiesPerCond <- function(e, expDesign){
	
	conditionMeanIntensities <- data.frame(row.names = rownames(e)) 
	
	skip = 1
	for(cond in names(expDesign)){
		
		replicates = expDesign[[cond]]
		#tag <- paste(cond,"MedianIntensity",sep="")
		tag = cond
		xCondFirstReplicate = skip
		xCondLastReplicate = xCondFirstReplicate+replicates-1
		xCondSpectralCount = data[,xCondFirstReplicate:xCondLastReplicate]
		meanIntensity = list()
		
		if(xCondFirstReplicate == xCondLastReplicate){
			meanIntensity[[tag]] = e[,xCondFirstReplicate]
		}else{
			meanIntensity[[tag]] = apply(e[,xCondFirstReplicate:xCondLastReplicate], 1 ,FUN=mean)
		}
		conditionMeanIntensities = data.frame(conditionMeanIntensities,meanIntensity)
		skip = skip+replicates
		
	}
	return(conditionMeanIntensities)
}

calculateRatios  <- function(conditionMedianIntensities){
	
	condNames = names(conditionMedianIntensities)
	medianRatios = data.frame(row.names = rownames(conditionMedianIntensities)) 
	for(i in 2:length(condNames)){
		
		ratio = list()
		ratio[[condNames[i]]] = conditionMedianIntensities[,i]/conditionMedianIntensities[,1]	
		medianRatios = data.frame(medianRatios,ratio)
	}
	
	return(medianRatios)
} 


createExpressionDataset <- function(logNormData,expDesign){
	### create expressionDataset
	
	conditions = c()
	
	for(cond in 1:length(names(expDesign))){
		conditions = c(conditions,rep(cond,expDesign[[as.character(cond)]]))
	}
	
	pData = data.frame(conditions)
	names(pData) = c("condition")
	
	#print(pData)
	#print(names(pData))
	
	rownames(pData) = colnames(logNormData)
	
	
	metadata = data.frame(labelDescription = "conditon label",row.names= names(pData))
	phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata) 
	annotation = "test"
	
	experimentData = new("MIAME", name = "Erik Ahrne",
			lab = "Proteomics Core Facillity, Biozentrum, Basel", contact = "erik.ahrne@unibas.ch",
			title = "test proteomics bioC", abstract = "An example ExpressionSet",
			url = "www.lab.not.exist", other = list(notes = "Created from text files"))
	
	eset = new("ExpressionSet", exprs = logNormData,
			phenoData = phenoData, experimentData = experimentData,
			annotation = annotation)
	
	### create expressionDataset end
}

spectralCountsPerCondition <- function(data,expDesign,spectralCountColumnsStart, proteins){
	
	spectralCounts = data.frame(row.names = proteins) 
	
	skip =0
	for(cond in names(expDesign)){
		
		replicates = expDesign[[cond]]
		tag = paste(cond,"SpectralCount",sep="")
		xCondFirstReplicate = spectralCountColumnsStart+skip
		xCondLastReplicate = xCondFirstReplicate+replicates-1
		xCondSpectralCount = data[,xCondFirstReplicate:xCondLastReplicate]
		spectralCountSum = list()
		
		if(xCondFirstReplicate == xCondLastReplicate){
			spectralCountSum[[tag]] = data[,xCondFirstReplicate]
		}else{
			spectralCountSum[[tag]] = apply(data[,xCondFirstReplicate:xCondLastReplicate], 1 ,FUN=sum)
		}
		spectralCounts = data.frame(spectralCounts,spectralCountSum)
		skip = skip+replicates
	}
	
	return(spectralCounts)
}

hClustHeatMap <- function(eset, condNames= condNames, orgCondNames= orgCondNames ,main=""){
	
	data = exp(exprs(eset))
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
	
	### pimped 1
	#color.map = function(cond) {if(cond==condNames[1]) "#FF0000" else "#0000FF" }
	color.map = function(cond) {colors()[10*which(cond== condNames)] }
	samplecolors = unlist(lapply(eset$condition,color.map))
	
	
	heatmap.2(normPerProteinData
			, col=topo.colors(75)
			, scale="none"
			, key=TRUE
			, symkey=FALSE
			, density.info="none"
			, trace="none"
			, cexRow=0.5
			, cexCol = 0.9
			,ColSideColors=samplecolors
			,Rowv=as.dendrogram(proteins.tree)
			,Colv=as.dendrogram(msrun.tree)
			,main=main
	)
	
	
	#legend("topright", c("control","other"), fill=c("#FF0000","#0000FF"))
	legend("topright", orgCondNames, fill=colors()[10*c(1:length(condNames))])
}

getQvals <- function(scores, isDecoy){
	
	
	targetScores = scores[!isDecoy]
	decoyScores = scores[isDecoy]
	
	qvals = c()
	for(score in scores){
		
		qval = sum(decoyScores >= score) / sum(targetScores >= score)
		qvals = c(qvals,qval)
	}
	
	return(qvals)
	
}

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

pairs.annot<-
		function(normData, main="", diagText=c(),...) {
	
	
	count <- 1
	
	panel.lm <-
			function (x, y, col = par("col"), bg = NA, pch = par("pch"),
					cex = 1, col.lm = "red", lwd=par("lwd"), ...)
	{
		points(x, y, pch = pch, col = col, bg = bg, cex = cex)
		
		ok = is.finite(x) & is.finite(y)
		if (any(ok))
			abline(lm(y~x,subset=ok), col = col.lm, ...)
	}
	
	panel.sse <-
			function(y, x, digits=2)
	{
		
		usr = par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		
		model = summary(lm(y~x))
		r2= model$r.squared
		#r=sqrt(r2)*sign(model$coef[2,1])
		#p= model$coef[2,4]
		
		txt = round(r2, digits)
		txt = bquote(R^2*"" == .(txt))
		
		text(0.5, 0.5, txt, cex=1.8)
		
		#txt = round(r2, digits)
		#txt = bquote(r^2 == .(txt))
		#text(0.5, 0.5, txt, cex=1.5)
		
		#txt = round(p, digits)
		#txt = bquote(P == .(txt))
		#text(0.5, 0.3, txt, cex=1.5)
		
	}
	
	panel.txt <- function(x, y, labels, cex, font,digits=2, ...){
		
		txt = colnames(normData)[count]
		#txt = bquote(r == .(txt))
		text(0.5, 0.6, txt, cex=1.8)
		
		if(length(diagText) > 0){
			txt = round(diagText[count], digits)
			txt = bquote(gf == .(txt))
			text(0.5, 0.4, txt, cex=2)
		}
		count <<-count+1
	}
	
	pairs(normData,lower.panel=panel.sse,upper.panel=panel.lm, text.panel=panel.txt,main=main, ...)
}

isCon <- function(accessions){
	
	#a <- c( "ACON_MOUSE","CON_MOUSE",":CON_MOUSE", "PROTEIN_CON" )
	#regexpr("\\Wcon_",a,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("^con_",a,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("_con",a,ignore.case=TRUE, perl=TRUE) > -1 
		
	return(regexpr("\\Wcon_",accessions,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("^con_",accessions,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("_con",accessions,ignore.case=TRUE, perl=TRUE) > -1 )
	
}

isDecoy <-function(accessions){
	
	return(((regexpr("rev_",accessions,ignore.case=TRUE) > -1)
						| (regexpr("decoy_",accessions,ignore.case=TRUE) > -1)
						| (regexpr("reverse_",accessions,ignore.case=TRUE) > -1)))
	
}

plotDeFDR <- function(fit,condition.names ,qval.cutoffs=seq(0,1,0.01) ){
	
	sign.prot = data.frame(row.names=qval.cutoffs)
	e.fp = data.frame(row.names=qval.cutoffs)
	
	i = 2
	for(cn in condition.names[2:length(condition.names)]){
		
		qvalues = p.adjust(fit$p.value[,i],method="fdr")
		sign.DE.proteins = c()
		expected.fp = c()
		for(fdr in qval.cutoffs){
			
			sign = sum(qvalues <= fdr )
			sign.DE.proteins = c(sign.DE.proteins,sign)
			expected.fp = c(expected.fp,sign*fdr)
		}
		
		s.tmp = list()
		e.tmp = list()
		
		s.tmp[[cn]] = sign.DE.proteins
		e.tmp[[cn]] = expected.fp
		
		sign.prot = data.frame(sign.prot,s.tmp)
		e.fp = data.frame(e.fp,e.tmp)
		
		i = i +1 
	}
	
	maxDispFdr = 0.1
	ylimExpectedFp = c(min(e.fp[qval.cutoffs < maxDispFdr,])-1,max(e.fp[qval.cutoffs < maxDispFdr,])+1) 
	limSignProt = c(min(sign.prot[qval.cutoffs < maxDispFdr,])-1,max(sign.prot[qval.cutoffs < maxDispFdr,])+1) 
	xlim = c(0,0.11)
	
#	par(mfrow=c(2,1))
	
	plot(0,0,ylim = limSignProt, xlim=xlim, type="n", xlab="False Discovery Rate", ylab="Sign. Features Counts")
	i = 1
	for(cn in condition.names[2:length(condition.names)]){
		lines(qval.cutoffs,sign.prot[,i],col=i)
		i = i + 1
	}
	legend("bottomright",condition.names[2:length(condition.names)],fill=1:(length(condition.names)-1))
	
#	
#	plot(0,0,ylim = ylimExpectedFp, xlim=limSignProt, type="n", xlab="Sign. Features Counts", ylab="Expected False Positive Counts")
#	i = 1
#	for(cn in condition.names[2:length(condition.names)]){
#		lines(sign.prot[,i],e.fp[,i],col=i)
#		i = i + 1
#	}
#	legend("topleft",condition.names[2:length(condition.names)],fill=1:(length(condition.names)-1))
#	
#	
	#par(mfrow=c(1,1))
	
}

plotValidPvalCutOff <- function(fit,condition.names ,qval.cutoffs=seq(0,1,0.01) ){
	
	sign.prot = data.frame(row.names=qval.cutoffs)
	e.fp = data.frame(row.names=qval.cutoffs)
	
	i = 2
	for(cn in condition.names[2:length(condition.names)]){
		
		pvalues = fit$p.value[,i]
		sign.DE.proteins = c()
		expected.fp = c()
		for(fdr in qval.cutoffs){
			
			sign = sum(pvalues <= fdr )
			sign.DE.proteins = c(sign.DE.proteins,sign)
			expected.fp = c(expected.fp,sign*fdr)
		}
		
		s.tmp = list()
		e.tmp = list()
		
		s.tmp[[cn]] = sign.DE.proteins
		e.tmp[[cn]] = expected.fp
		
		sign.prot = data.frame(sign.prot,s.tmp)
		e.fp = data.frame(e.fp,e.tmp)
		
		i = i +1 
	}
	
	maxDispFdr = 0.1
	ylimExpectedFp = c(min(e.fp[qval.cutoffs < maxDispFdr,])-1,max(e.fp[qval.cutoffs < maxDispFdr,])+1) 
	limSignProt = c(min(sign.prot[qval.cutoffs < maxDispFdr,])-1,max(sign.prot[qval.cutoffs < maxDispFdr,])+1) 
	xlim = c(0,0.11)
	
#	par(mfrow=c(2,1))
	
	plot(0,0,ylim = limSignProt, xlim=xlim, type="n", xlab="False Discovery Rate", ylab="Sign. Features Counts")
	i = 1
	for(cn in condition.names[2:length(condition.names)]){
		lines(qval.cutoffs,sign.prot[,i],col=i)
		i = i + 1
	}
	legend("bottomright",condition.names[2:length(condition.names)],fill=1:(length(condition.names)-1))
	
#	
#	plot(0,0,ylim = ylimExpectedFp, xlim=limSignProt, type="n", xlab="Sign. Features Counts", ylab="Expected False Positive Counts")
#	i = 1
#	for(cn in condition.names[2:length(condition.names)]){
#		lines(sign.prot[,i],e.fp[,i],col=i)
#		i = i + 1
#	}
#	legend("topleft",condition.names[2:length(condition.names)],fill=1:(length(condition.names)-1))
#	
#	
	#par(mfrow=c(1,1))
	
}


plotVolcanoes <- function(fit, medianRatios,condition.names,conditionNas,conditionLargeSd,ratioCutOff=ratioCutOff, qvalVolcano=qvalVolcano, deFdrCutoff=deFdrCutoff,BCutoff=BCutoff){
	
	controlNA = conditionNas[,1]
	controlLargeSd = conditionLargeSd[,1]
	
	xlim <- c(min(log(medianRatios))*0.95, max(log(medianRatios))*1.05)
	
	for(i in 2:length(condition.names)){
		
		condNA = conditionNas[,i]
		condLargeSd = conditionLargeSd[,i]
		nas = controlNA | condNA
		number = 15
		
		tt = topTable(fit, number=length(fit$p.value[,1]), coef=i, adjust.method="fdr")
		
		order = as.integer(rownames(tt))
		
		
		
#		 ### mean FC
#		if(meanFC){
#			plot(tt$logFC,tt$B, pch=20, xlab="Log Fold Change", ylab="Log Odds",col="grey", main=condition.names[i])
#			text(tt$logFC[indexTop],tt$B[indexTop],tt$ID[indexTop],pos=3,col="grey", cex=0.5)
#			points(fit$coefficients[,i][indexValid],fit$lods[,i][indexValid],pch=20 ,col="black")
#			points(fit$coefficients[,i][condLargeSd],fit$lods[,i][condLargeSd] ,pch=2,col="red")
#			points(fit$coefficients[,i][controlLargeSd],fit$lods[,i][controlLargeSd] ,pch=6,col="red")
#			points(fit$coefficients[,i][condNA],fit$lods[,i][condNA],pch=2 ,col="blue")
#			points(fit$coefficients[,i][controlNA],fit$lods[,i][controlNA],pch=6 ,col="blue")
#			labels = c("Valid","Control NAs","Cond. NAs","Control large SD","Cond. large SD")
#			legend("bottomleft",labels, col=c("black", "blue","blue","red","red"), pch= c(20,2,6,2,6))
#			
#		}else{
		### median FC
						
		### qvals
		if(qvalVolcano){
			qvals = abs(log10(p.adjust(fit$p.value[,i],method="fdr")))
			
			#plot(sort(qvals),sort(abs(log10(tt$adj.P.Val))))
			
			validCond = (tt$adj.P.Val < deFdrCutoff) & abs(log(medianRatios[,i-1])[order]) > log(ratioCutOff)  
			#validCond = (tt$adj.P.Val < deFdrCutoff) & (abs(tt$logFC) > log(ratioCutOff))  
			indexValid = as.integer(rownames(tt)[validCond])
			indexTop =  order(tt$adj.P.Val[validCond], decreasing=FALSE)[1:number]
			
			plot(log(medianRatios[,i-1])[order],abs(log10(tt$adj.P.Val)), pch=20, xlab="Log Fold Change", ylab="abs log10 qval",col="grey", main=paste("Median Ratios",condition.names[i]), xlim=xlim)
			text(log(medianRatios[,i-1][order][indexTop]),abs(log10(tt$adj.P.Val[indexTop])),tt$ID[indexTop],pos=3,col="grey", cex=0.5)
			points(log(medianRatios[,i-1][indexValid]),qvals[indexValid],pch=20 ,col="black")
			points(log(medianRatios[,i-1][condLargeSd]),qvals[condLargeSd] ,pch=2,col="red")
			points(log(medianRatios[,i-1][controlLargeSd]),qvals[controlLargeSd] ,pch=6,col="red")
			points(log(medianRatios[,i-1][condNA]),qvals[condNA],pch=2 ,col="blue")
			points(log(medianRatios[,i-1][controlNA]),qvals[controlNA],pch=6 ,col="blue")
			labels = c("Valid","Control NAs","Cond. NAs","Control large SD","Cond. large SD")
			legend("bottomleft",labels, col=c("black", "blue","blue","red","red"), pch= c(20,2,6,2,6))
		}else{
			
			### lods
			
			validCond = (tt$B > BCutoff)  & abs(log(medianRatios[,i-1])[order]) > log(ratioCutOff) 
			validCond = (tt$B > BCutoff) & (abs(tt$logFC) > log(ratioCutOff))  
			indexValid = as.integer(rownames(tt)[validCond])
			indexTop =  order(tt$B[validCond], decreasing=TRUE)[1:number]
			
			plot(log(medianRatios[,i-1])[order],tt$B, pch=20, xlab="Log Fold Change", ylab="Log Odds",col="grey", main=paste("Median Ratios",condition.names[i]), xlim=xlim)
			text(log(medianRatios[,i-1][order][indexTop]),tt$B[indexTop],tt$ID[indexTop],pos=3,col="grey", cex=0.5)
			points(log(medianRatios[,i-1][indexValid]),fit$lods[,i][indexValid],pch=20 ,col="black")
			points(log(medianRatios[,i-1][condLargeSd]),fit$lods[,i][condLargeSd] ,pch=2,col="red")
			points(log(medianRatios[,i-1][controlLargeSd]),fit$lods[,i][controlLargeSd] ,pch=6,col="red")
			points(log(medianRatios[,i-1][condNA]),fit$lods[,i][condNA],pch=2 ,col="blue")
			points(log(medianRatios[,i-1][controlNA]),fit$lods[,i][controlNA],pch=6 ,col="blue")
			
			prob = round(exp(BCutoff)/(1+exp(BCutoff)),2)
			labels = c(paste("Prob. DE >", prob ),"Control NAs","Cond. NAs","Control large SD","Cond. large SD")
			legend("bottomleft",labels, col=c("black", "blue","blue","red","red"), pch= c(20,2,6,2,6))
			
		}											
	}
	
#	if(opt$verbose){
#		if(meanFC){
#			print("Mean based ratios displayed in volcano plot") 
#		}else{
#			print("Median based ratios displayed in volcano plot") 
#		}
#	}
	
}


### FUNCTIONS
