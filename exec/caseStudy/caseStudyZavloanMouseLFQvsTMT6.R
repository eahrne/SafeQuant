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


# plot uncorrected protein ratio correlation
plotXYDensity(sqaLFQProtein$ratio[,1], sqaTMTProtein$ratio[ gsub("_$","",rownames(sqaLFQProtein$ratio)) ,],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

# plot uncorrected spectrum  ratios vs. reference ratios
plotXYDensity(sqaLFQProtein$ratio[ paste(fData(sqaTMTSpectrum$eset)$proteinName,"_",sep="") ,], sqaTMTSpectrum$ratio[,1],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )


#@TODO plotXYDensity(apply(sqaTMTSpectrum$), sqaLFQProtein$ratio[ paste(fData(sqaTMTSpectrum$eset)$proteinName,"_",sep="") ,] - sqaTMTSpectrum$ratio[,1],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )


head(esetLFQProtein)

esetLFQPeptide

esetLFQSpectrum


unique(fData(esetLFQPeptide)$proteinName)

gsub("_","",rownames(esetLFQPeptide)) %in% rownames(esetTMTPeptide)


unique(fData(esetLFQPeptide)$proteinName)



intersect(gsub("_$","",rownames(esetLFQProtein)), rownames(esetTMTProtein))





sqaLFQProtein$ratio

sqaTMTProtein$ratio


