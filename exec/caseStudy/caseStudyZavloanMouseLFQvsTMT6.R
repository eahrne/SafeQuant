# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

source("/Users/erikahrne/dev/R/workspace/SafeQuant/exec/caseStudy/initCaseStudySession.R")


### LOAD DATA

load("/Users/erikahrne/dev/R/workspace/SafeQuant/data/zavolanMouseLFQvsTMT6.rda")
ls(pattern="*eset")
### same protein annotations for TMT and LFQ
#fData(esetLFQSpectrum)$proteinName <- fData(esetTMTPeptide[as.character(fData(esetLFQSpectrum)$peptide),])$proteinName

### LOAD DATA END

# re-run createZavolanMouseLFQTvsTMT6estData.R (2 peptide per protein cutoff)

# get protein level reference ratios from LFQ data

# create protein sqa objects
sqaLFQProtein <- safeQuantAnalysis(esetLFQProtein,method="global")
sqaTMTProtein <- safeQuantAnalysis(esetTMTProtein,method="global")

# create spectrum sqa objects
sqaTMTSpectrum <- safeQuantAnalysis(esetTMTSpectrum,method="global")

ratioTMT <-  sqaTMTSpectrum$ratio[,1]


### NO FILTERING
ratioLFQ <- sqaLFQProtein$ratio[ as.character(fData(sqaTMTSpectrum$eset)$proteinName) ,]

# calculate interferemce level , Rref = LFQ2/LFQ1 , il = (TMT2 - Rref*TMT1) / (1-Rref) 
tmtMedianPerCond <-  getSignalPerCondition(sqaTMTSpectrum$eset)
tmtMinInt <- apply(sqaTMTSpectrum$eset,1,min)
tmtInterference <- (tmtMedianPerCond[,2]- (2^ratioLFQ)*tmtMedianPerCond[,1]) / (1-(2^ratioLFQ))

tmpSel <- is.na(tmtInterference)
tmtInterference[tmpSel] <- 0
tmtInterference[(tmtInterference > tmtMinInt) |  (tmtInterference < 0)] <- tmtMinInt[(tmtInterference > tmtMinInt) |  (tmtInterference < 0)  ]
tmtInterference[tmpSel] <- NA

#ratioDelta <- (abs(ratioLFQ) - abs(ratioTMT))*sign(ratioLFQ)*sign(ratioTMT)

## plot minInt vs. interference
plotXYDensity(log2(tmtMinInt)[ (tmtInterference != tmtMinInt) & (abs(sqaTMTSpectrum$ratio[,1]) > 2)], log2(tmtInterference)[ (tmtInterference != tmtMinInt) & (abs(sqaTMTSpectrum$ratio[,1]) > 2)])
#abline(h=0, lty=2)

plotXYDensity(log2(tmtMinInt)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 0)], log2(tmtInterference)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 0)]) 
# gloally 55% of tmtMinInt is noise

plotXYDensity(log2(tmtMinInt)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 1)], log2(tmtInterference)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 1)])
plotXYDensity(log2(tmtMinInt)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 2)], log2(tmtInterference)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 2)])
plotXYDensity(log2(tmtMinInt)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 2.5)], log2(tmtInterference)[ (tmtInterference != tmtMinInt) & !is.na(ratioLFQ) & (abs(ratioLFQ) > 2.5)])

print("DONE")

plot(c(1:10),(c(1:10))/0.3)

plot(c(1:10)-1,(c(1:10) - 1)/0.3)



tmtMedianPerCondAdjusted <- tmtMedianPerCond - (tmtMinInt * 0.7) 


ratioTMTAdjusted <- log2(tmtMedianPerCondAdjusted[,2] / tmtMedianPerCondAdjusted[,1])

plotXYDensity(ratioTMT,ratioTMTAdjusted)
abline(coef=c(0,1))


### HERE
plotXYDensity(ratioLFQ,ratioTMT)
plotXYDensity(ratioLFQ,ratioTMTAdjusted)

### DIAGNOSTIC PLOTS

