# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

source("/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/initCaseStudySession.R")


### LOAD DATA

load("/Users/erikahrne/dev/R/workspace/SafeQuant/data/zavolanMouseLFQvsTMT6.rda")
ls(pattern="*eset")
### same protein annotations for TMT and LFQ
fData(esetLFQSpectrum)$proteinName <- fData(esetTMTPeptide[as.character(fData(esetLFQSpectrum)$peptide),])$proteinName

### LOAD DATA END

# re-run createZavolanMouseLFQTvsTMT6estData.R (2 peptide per protein cutoff)

# get protein level reference ratios from LFQ data

# create protein sqa objects
sqaLFQProtein <- safeQuantAnalysis(esetLFQProtein,method="global")
sqaTMTProtein <- safeQuantAnalysis(esetTMTProtein,method="global")

# create spectrum sqa objects
sqaTMTSpectrum <- safeQuantAnalysis(esetTMTSpectrum,method="global")

# @TODO select for proteins significantly regulated in LFQ

# plot uncorrected protein ratio correlation
plotXYDensity(sqaLFQProtein$ratio[,1], sqaTMTProtein$ratio[ gsub("_$","",rownames(sqaLFQProtein$ratio)) ,],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

# plot uncorrected spectrum  ratios vs. reference ratios
plotXYDensity(sqaLFQProtein$ratio[ paste(fData(sqaTMTSpectrum$eset)$proteinName,"_",sep="") ,], sqaTMTSpectrum$ratio[,1],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

# plot log2 ratio difference vs tmt intensity
plotXYDensity(log2(apply(exprs(sqaTMTSpectrum$eset),1,min))
	, sqaLFQProtein$ratio[ paste(fData(sqaTMTSpectrum$eset)$proteinName,"_",sep="") ,] - sqaTMTSpectrum$ratio[,1]
	, xlab="log2 TMT intensity sum"
	, ylab="(LFQ log2 Ratio) - (TMT log2 Ratio) " )
abline(h=0, lty=2)


boxplot(sqaLFQProtein$ratio[ paste(fData(sqaTMTSpectrum$eset)$proteinName,"_",sep="") ,] - sqaTMTSpectrum$ratio[,1])
abline(h=0, lty=2)

plotVolcano(sqaLFQProtein, adjusted=F)
