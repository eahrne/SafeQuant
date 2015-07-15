# TODO: Add comment
# 
#
# If you want to disable tooltips in Firefox, start by typing about:config into the address bar.
# In the list that appears, find browser.chrome.toolbar_tips. An easy way to find it is to start typing browser.chro into the filter box.
# Once you have found the item, double-click on it to change it from true to false. Now when you are browsing,
# tooltips should be disabled. Simply double-click again to change things back.
### Author: erikahrne
###############################################################################

### call volcanoPlotSVG for all vs. control condition pairs
plotAllVolcanoesSVG <- function(ratiosMedianNormIntPerCond,qvaluesPerCond,cvsPerCond
		,controlCondition=controlCondition
		,nonControlConditions=nonControlConditions
		,qvalueCutOff =0.01
		,ratioCutOff=2
		,outPath
		,labels1=rownames(qvaluesPerCond)
		,labels2=rep("",nrow(qvaluesPerCond))
		, ...){
	
	fList <- list()
	### iterate over all nonControlConditions
	c <- 1
	for(cond in nonControlConditions){

		qvalues <- qvaluesPerCond[,cond]
		ratios <- ratiosMedianNormIntPerCond[,cond]
		cvMaxPerPair <-  apply(cvsPerCond[, names(cvsPerCond)  %in% c(cond,controlCondition) ],1,FUN=function(cvs){return(max(cvs, na.rm=T))})
		
		cond <- gsub("[ 	]","",cond)
		#svgFileName <- paste(cond,".svg",sep="")
		svgFileName <- paste("volcanoplot",c,".svg",sep="") ## name svg file volcanoplot1.svg etc. 
		svgFileNamePath <- paste(outPath,svgFileName,sep="")
		
		volcanoPlotSVG(ratios,qvalues,cvMaxPerPair,condName=cond ,controlCondName=controlCondition, qvalueCutOff = qvalueCutOff,ratioCutOff=ratioCutOff, labels1=labels1,labels2=labels2 ,svgFileName=svgFileNamePath, ... )
		
		### path relative to .html file
		fList[svgFileName] <- paste("volcano plot ",cond," vs ",controlCondition,sep="")
		c <- c+1
	}
	
	return(fList)
	
}

volcanoPlotSVG <- function(ratios,qvalues,cvMaxPerPair,svgFileName=svgFileName,condName=condName,labels1=rep("",length(ratios)),labels2=rep("",length(ratios)) ,controlCondName=controlCondName,qvalueCutOff=0.05,ratioCutOff=2,higlightSel=rep(0,length(ratios)), ...  ){
	
	devSVGTips(svgFileName,toolTipMode=2)
	cvColorPalette <- rev(rich.colors(32))
	
	### plot extras params
	ratioCutOffAbsLog2 <- abs(log2(ratioCutOff))
	CVcolors <- cvColorPalette[round(1+ cvMaxPerPair/max(cvMaxPerPair) *31)] 
	
	qvaluesAbsLog <- abs(log10(qvalues)) 
	ratiosLog <- log2(ratios)
	
	par(fig=c(0,1,0.18,1), mar=c(5.1,2.1,4.1,15.1))
	plot(ratiosLog,qvaluesAbsLog,xlab= paste("log2(",condName,"/",controlCondName,")", sep="" ), type="n", ...)
	mtext(side=4,line=1,"log10 p-value")
	
	
	
	### add tooltip
	invisible(sapply(1:length(ratiosLog), function(i)
					{	setSVGShapeToolTip(title=labels1[i]
								,desc1=labels2[i]
								,desc2= paste("ratio:",signif(ratios[i],2),"p-value:",signif(qvalues[i],2),sep=" ")
						
						)
						setSVGShapeURL(paste("http://www.google.com/search?q=",labels1[i],"&btnI",sep="")) ### guess link	
    					points(ratiosLog[i], qvaluesAbsLog[i], pch=20, col=CVcolors[i])}))
	
	### draw valid squares contianing valid proteins/peptides 
	lines(c(-ratioCutOffAbsLog2,-ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),1000), col="grey")
	lines(c(ratioCutOffAbsLog2,ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),1000), col="grey")
	lines(c(-1000,-ratioCutOffAbsLog2),c(abs(log10(qvalueCutOff)),abs(log10(qvalueCutOff))), col="grey")
	lines(c(ratioCutOffAbsLog2,1000),c(abs(log10(qvalueCutOff)),abs(log10(qvalueCutOff))), col="grey")
	
	### highlight selection w/ black circle
	points(ratiosLog[higlightSel],qvaluesAbsLog[higlightSel])
	
	### add heat color bar
	par(fig=c(0,1,0,0.25), new=TRUE)
	cvColorstrip(cvColorPalette, maxCV=max(cvMaxPerPair)*100)
	par(mfrow=c(1,1))
	
	dev.off()
	
}


### CREATE SVG IMAGE OF hClustHeatMap3 plot
heatMap2SVG <-  function(dm,addLink=T, svgFileName=svgFileName, isUniprot=NA,desc=rep("NA",nrow(dm)), ...){
	
	### @TODO if isUniprot and addLink -> link to uniprot 	
	### if addLink and not uniprot -> use google to guess
	
	devSVGTips(svgFileName,toolTipMode=2, toolTipFontSize = 8,width = 12, height = 12)
	
	#hm <- heatmap.2(dm,trace="none", ...)
	
	hm <- hClustHeatMap3(dm
			, margins = c(20,20) ### keep tool tips within plot marginfs (so they won't )
			,keysize=1.1
			,...
					
					
					
	)
	
	rowLabels <- 1:nrow(dm)
	
	if(!is.na(colnames(hm$carpet)[1])){
		rowLabels <- colnames(hm$carpet)
	}
	
	desc2 <- apply(t(hm$carpet),1,FUN=function(t){ return(paste(round(t),collapse=","))   }) ### string showing up/down reg e.g -2,1,0
	drawHeatMapToolTips(numRow = nrow(dm),addLink=addLink,rowLabels=rowLabels,desc1=desc,desc2=desc2)
	
	dev.off()
}

drawHeatMapToolTips <- function(numRow = numRow,addLink=addLink,rowLabels=rowLabels,desc1=desc1,desc2=desc2, margins=1){
	### +c(9.9,-2.75,-2.6,+9.9) there to afjust for margins = c(20,20) & keysize 1.1
	if(margins == 0 ){
		#par(mar=c(4.1,12.8,12.8,4.1)+c(9.9,+4.45,+5.5,+9.9)) ### adjust plot margins of "invisible plot" to match hm2 plotting area
		par(mar=c(4.1,21.8,21.8,4.1)) ### adjust plot margins of "invisible plot" to match hm2 plotting area
	}else if(margins == 1){
		#par(mar=c(3.3,12.8,14,3.3)) ### adjust polt margins of "invisible plot" to match hm2 plotting area (WITH COLSIDE COLORS)
		par(mar=c(3.3,12.8,14,3.3)+c(9.9,+4.45,+5.6,+9.9)) ### adjust plot margins of "invisible plot" to match hm2 plotting area (WITH COLSIDE COLORS)
	}else if(margins == 2){
		par(mar=c(3.3,13.9,14,3.3)+c(9.9,+4.45,+5.5,+9.9)) ### adjust plot margins of "invisible plot" to match hm2 plotting area (WITH ROWSIDE COLORS)
	}
	
	seqrchQuery <- gsub("^[tr sp]{,2}\\|","",rowLabels)
	dim <- 100
	plot(0:dim,0:dim, type="n",xaxt="n",yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
	
	px <- dim/numRow
	indent <- 10^-10		
	
	for(i in 0:(numRow-1)){
		
		#setSVGShapeToolTip(title=paste(rowLabels[i+1]),desc=desc[i+1])
		setSVGShapeToolTip(title=paste(rowLabels[i+1])
				,desc1=desc1[i+1]
				,desc2=desc2[i+1]
			)
		#print(desc[i+1])
		if(addLink){
			setSVGShapeURL(paste("http://www.google.com/search?q=",seqrchQuery[i+1],"&btnI",sep="")) ### guess link
		}
		
		ybottom <- px*i
		ytop <- (px*i) + px
		
		
		if(ybottom == 0){
			ybottom  <- indent
		}
		if(ytop == dim){
			ytop  <- ytop - indent
		}
		
		#print(paste(i,ybottom))
		
		rect(indent,ybottom,dim-indent ,ytop,border = "transparent")
		#rect(indent,ybottom,10-indent ,ytop)
	}
	
}

svgLinkHTML <- function(file,desc){
	return(paste('<p><a href="',file,'">',desc,'</a></p>', sep=""))
}

svgEmbedHTML <- function(file,desc){
	return(paste('<object data="',file,'" type="image/svg+xml">',desc,'</object>', sep=""))
}

### create simple HTML with links to svg file listed in fileList
createSVGLinksHTML <- function(fileList, htmlFile, isEmbedded=F){
	
	headerText <-'<!DOCTYPE html>
			<html lang="en">
			<head>
			<meta charset="UTF-8" />
			<title>SVG PLOTS</title>
			</head>
			<body>
			<center>
			'
	footerText 	<- '
			</center>
			</body>
			</html>'
	
	linksHTML <- c()
	for(file in names(fileList)){
		if(isEmbedded){
			linksHTML <- c(linksHTML,svgEmbedHTML(file,fileList[file]))
		}else{
			linksHTML <- c(linksHTML,svgLinkHTML(file,fileList[file]))
		}
	}
	
	### print to html file
	cat(headerText
			,linksHTML
			,footerText
			, file=htmlFile )
	
}

plotXYSVG <- function(x,y,thumbLabel=rep("",length(x)), thumbDesc1=rep("",length(x)), fileName=fileName, ...){
	
	### discard NAs
	ok <- is.finite(x) & is.finite(y)
	thumbLabel <- thumbLabel[ok]
	thumbDesc1 <- thumbDesc1[ok]
	x <- x[ok]
	y <- y[ok]
	
	#devSVGTips(fileName,toolTipMode=2)
	devSVGTips(fileName,toolTipMode=2, width=10, height=10)
	
	### leave a bit of margin so that thumbnail doesn't overflow
	#par(fig=c(0,1,0.18,1), mar=c(5.1,2.1,4.1,15.1))
	par(fig=c(0,1,0.18,1), mar=c(5.1,4.1,4.1,15.1))
	plot(x,y, type="n", ...)
	#mtext(side=4,line=1,"bla")
	
	### add tooltip
	invisible(sapply(1:length(x), function(i)
					{	
						setSVGShapeToolTip(
								title=thumbLabel[i]
								,desc1=thumbDesc1[i]
								
						)
						setSVGShapeURL(paste("http://www.google.com/search?q=",thumbDesc1[i],"&btnI",sep="")) ### guess link	
    					points(x[i], y[i], pch=19)}))
	
}