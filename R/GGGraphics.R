
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
	
	# thresholds		
	#	 pvalue thrs
	p =  p + geom_abline(intercept = -log10(pValueThrs),slope=0, lty=thrsLineLty, col=thrsLineCol)
	# ratio thrs
	p =  p + geom_vline(xintercept=c(-log2RatioThrs,log2RatioThrs), lty=thrsLineLty, col=thrsLineCol)	   
	
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
	dfLab = subset(data, (pValue < pValueThrs) & (abs(ratio) > log2RatioThrs))
	if(nrow(dfLab) > 0){
		p = p + geom_text_repel(data= dfLab,aes(x =ratio,y=-log10(pValue),label=geneName ))
	}
	
	return(p)
	
}