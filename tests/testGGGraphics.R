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

cv = apply(sqa$cv[,1:2],1,max, na.rm=T)*100
ggDf = data.frame(ratio = sqa$ratio[,1], pValue = sqa$pValue[,1], geneName = fData(sqa$eset)$proteinName ,  ac=NA, cv = cv, description=fData(sqa$eset)$proteinDescription  )



tmpFile <- paste(tempdir(),"/tmp_ggplot.pdf",collapse="",sep="")

pdf(tmpFile)


plot(ggVolcanoPlot(data=ggDf,pValueThrs=pValueThrs,log2RatioThrs=log2RatioThrs,thrsLineCol = "lightgrey",thrsLineLty = 2, title="test title"))



dev.off()

cat("CREATED FILE", tmpFile, "\n")


# p = ggVolcanoPlot(data=ggDf,pValueThrs=pValueThrs,log2RatioThrs=log2RatioThrs,thrsLineCol = "lightgrey",thrsLineLty = 2 )
# ggplotly(p)


# testGGVolcanoPlot <- function(){
#   
#   cat(" --- testGGVolcanoPlot --- \n")
#   stopifnot(T)
#   cat(" --- testGGVolcanoPlot: PASS ALL TEST --- \n")
#   
# }