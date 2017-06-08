
#' Plots volcano, data points colored by max cv of the 2 compared conditions
#' @param data data.frame
#' @param log2RatioThrs default log2(0.5)
#' @param pValueThrs default 0.01
#' @param thrsLineCol default "lightgrey"
#' @param defalut 2
#' @param xlab default "log2 ratio" 
#' @param ylab default "-log10 pValue"
#' @param title default no title
#' @param xlim xlim
#' @param ylim ylim
#' @param abline c("none","both","ratio","pvalue")
#' @params topNlabels default 10, label top proteins/peptides ordered by p-value
#' @param textSize default 20
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
ggVolcanoPlot = function(data=data
		, title=""
		, pValueThrs=0.05
		, log2RatioThrs=0.5849625
		, thrsLineCol = "lightgrey"
		, thrsLineLty = 2
		, xlab = "log2 ratio"
		, ylab = "-log10 pValue"
		, textSize = 20
		, xlim = range(data$ratio,na.rm=T)
		, ylim = range(-log10(data$pValue),na.rm=T)
		, abline = c("both")
		, topNlabels = 10
){
	
	# plotted data	
	p =  ggplot(data,aes(x = ratio,y=-log10(pValue)
					,label=pValue
					,label2=geneName
					,label3=ac
					,lable4=description
					,color=cv)) 
	
	# axis 
	p = p + labs(list(x=xlab, y=ylab, title=title))
	p = p + theme_bw()
	p = p + scale_x_continuous(limits =xlim)
	p = p + scale_y_continuous(limits =ylim)
	
	# abline		
	#	 pvalue thrs
	if(abline %in% c("pvalue","both") ) p =  p + geom_abline(intercept = -log10(pValueThrs),slope=0, lty=thrsLineLty, col=thrsLineCol)
	# ratio thrs
	if(abline %in% c("ratio","both") )p =  p + geom_vline(xintercept=c(-log2RatioThrs,log2RatioThrs), lty=thrsLineLty, col=thrsLineCol)	   
	
	# point style
	p = p + geom_point() 
	p = p + scale_colour_gradientn(colours=c("red","yellow","blue"), name="C.V. (%)" )
	
	# theme
	p = p + theme(text = element_text(size=textSize)
			#, axis.text.x = element_text(angle=0, hjust=1)
			, legend.position="right"
			, legend.direction="vertical"
			, legend.title = element_text(size=textSize*0.8)
			, legend.text=element_text(size=textSize*0.6)
	
	)
	
	# add labels
	# geom_GeomTextRepel() has yet to be implemented in plotly (status at v. 4.5.6).
	# disp geneName above thrs 
	labPvalueThrs = ifelse(topNlabels > 0,sort(data$pValue)[min(topNlabels,length(data$pValue))],0)
	dfLab = subset(data, (pValue <= min(labPvalueThrs,pValueThrs) ) & (abs(ratio) >= log2RatioThrs))
	dfLab = dfLab[1:min(10,nrow(dfLab)),]
	if(nrow(dfLab) > 0){
		p = p + geom_text_repel(data= dfLab,aes(x =ratio,y=-log10(pValue),label=geneName ))
	}
	
	return(p)
	
}


#' Plots volcano of all condition comparisons
#' @param sqa SafeQuantAnalysis object
#' @param isAdjusted (T/F) plot adjusted pvalues
#' @param see ggVolcanoPlot
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
plotAllGGVolcanoes = function(sqa, isAdjusted=T ,...){
  # plot all volcanoes
  ctrlCondition = pData(sqa$eset)$condition[pData(sqa$eset)$isControl][1] %>% as.character
  caseConditions = setdiff(pData(sqa$eset)$condition %>% unique, ctrlCondition)
  
  if(isAdjusted){
    allPValue = sqa$qValue
  }else{
    allPValue = sqa$pValue
  }
  
  xlim  = range(sqa$ratio, na.rm=T)
  ylim = range(abs(log10(allPValue)), na.rm=T)
  
  for(cond in caseConditions){
    
    # compile df
    cv = apply(sqa$cv[, c(ctrlCondition,cond) ],1,max, na.rm=T)*100
    ggDf = data.frame(ratio = sqa$ratio[,cond], pValue=allPValue[,cond] ,geneName = fData(sqa$eset)$geneName ,  ac=fData(sqa$eset)$ac, cv = cv, description=fData(sqa$eset)$proteinDescription  )
  
    #plot
    plot(ggVolcanoPlot(data=ggDf, xlab = paste("log2", cond,"/",ctrlCondition ),xlim=xlim,ylim=ylim,  ... ))
  }
}
