# Author: ahrnee-adm
###############################################################################

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
  setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END

pValueThrs = 0.01
log2RatioThrs=0.001
thrsLineCol = "lightgrey"
thrsLineLty = 2

intensity = apply(getSignalPerCondition(sqa$eset)[,1:2],1,max, na.rm=T)
ggDf = data.frame(ratio = sqa$ratio[,1], pValue = sqa$pValue[,1]
                  , geneName = fData(sqa$eset)$proteinName
                  , ac=NA
                  , intensity = intensity
                  , description=fData(sqa$eset)$proteinDescription
                  , naHighLightSel= c(rep(T,300),rep(F,600))
                  , naCat = factor(c(rep("ctrl",100), rep("case",100), rep("both",100), rep("none",600)), levels =  c("ctrl","case","both","none"))
                  )

#ggDf[ggDf$geneName == "REV_prot158",]

tmpFile <- paste(tempdir(),"/tmp_ggplot.pdf",collapse="",sep="")

pdf(tmpFile)

plot(ggVolcanoPlot(data=ggDf,pValueThrs=pValueThrs,log2RatioThrs=log2RatioThrs,thrsLineCol = "lightgrey",thrsLineLty = 2, title="test title"))

dev.off()

cat("CREATED FILE", tmpFile, "\n")


#plotAllGGVolcanoes(sqa)


