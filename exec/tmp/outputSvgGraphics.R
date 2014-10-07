# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

if(userOptions$isCreateSVGFigs){
	
	htmlFile <- paste(userOptions$outputDirPath,userOptions$resultsFileLabel,".html",sep="")
	
	if(exists("PROTEINDESCDIC")){ ###TMT
		proteinDescription <- PROTEINDESCDIC[rownames(UNGROUPEDNORMINTDATA),1]
	}else{ ### LFQ
		proteinDescription <- PROGOUTDATA$Description
	}
	
	### VOLCANO PLOTS
	if(userOptions$isVolcanoPlots){
		fList <- plotAllVolcanoesSVG(RATIOSMEDIANNORMINTPERCOND,QVALUESPERCOND,CVSPERCOND,outPath=userOptions$outputDirPath
			,controlCondition=CONTROLCONDITION
			,nonControlConditions=NONCONTROLCONDITIONS
			,labels1=rownames(QVALUESPERCOND)
			,labels2=proteinDescription
		
		)
	}
	
	### HIERARCHICAL CLUSTERING ratios
	if(userOptions$isHClustPlot ){
		
		###  No SVG heat map if more than 2500 features
		if(nrow(UNGROUPEDNORMINTDATA) < 2500 ){
			heatMap2SVG(UNGROUPEDNORMINTDATA,svgFileName=paste(userOptions$outputDirPath,"hClustHeatMap.svg",sep="")
					,expDesign=EXPDESIGNOBJ$expDesign
					,conditionColors=CONDITIONCOLORS 
					,controlIntensities=MEDIANNORMINTDATA[,CONTROLCONDITION]
					,desc=proteinDescription
			)
    		fList["hClustHeatMap.svg"] <- "HClust HeatMap"
		}
	}
	
	if(!is.na(fList)){
		createSVGLinksHTML(fList,htmlFile)
		cat("CREATED FILE ", htmlFile,"\n")
	}
	
}
