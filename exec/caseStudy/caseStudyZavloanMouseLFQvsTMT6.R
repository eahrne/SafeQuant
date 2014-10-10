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

ratioDelta <- (abs(ratioLFQ) - abs(ratioTMT))*sign(ratioLFQ)*sign(ratioTMT)

# plot uncorrected protein ratio correlation
plotXYDensity(sqaLFQProtein$ratio[,1], sqaTMTProtein$ratio[ gsub("_$","",rownames(sqaLFQProtein$ratio)) ,],xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

# plot uncorrected spectrum  ratios vs. reference ratios
plotXYDensity(ratioLFQ,ratioTMT,xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

## plot log2 ratio difference vs tmt intensity
plotXYDensity(log2(apply(exprs(sqaTMTSpectrum$eset),1,min))
	, ratioDelta
	, xlab="log2 TMT intensity sum"
	, ylab="Compression Factor" )
abline(h=0, lty=2)

boxplot(ratioDelta)
abline(h=0, lty=2)

# SELECT FOR PROTEINS SIGNIFICANTLY REGULATED IN LFQ


selProteins <- unique(c(rownames(sqaLFQProtein$pValue)[sqaLFQProtein$pValue[,1] <= 0.01],rownames(sqaTMTProtein$pValue)[sqaTMTProtein$qValue[,1] <= 0.01]))

selLFQ <- rownames(sqaLFQProtein$pValue) %in% selProteins
#sqaLFQProteinSel <- .filterSQA(sqaLFQProtein,selLFQ)


ratioLFQSEL <- .filterSQA(sqaLFQProtein,selLFQ)$ratio[ as.character(fData(sqaTMTSpectrum$eset)$proteinName) ,]

ratioDeltaSel <- (abs(ratioLFQSEL) - abs(ratioTMT))*sign(ratioLFQSEL)*sign(ratioTMT)

# plot uncorrected spectrum  ratios vs. reference ratios
plotXYDensity(ratioLFQSEL,ratioTMT,xlab="LFQ log2 Ratio", ylab="TMT log2 Ratio" )

# plot log2 ratio difference vs tmt intensity
plotXYDensity(log2(apply(exprs(sqaTMTSpectrum$eset),1,sum))
		, ratioDeltaSel
		, xlab="log2 TMT intensity sum"
		, ylab="Compression Factor" )
abline(h=0, lty=2)

# plot log2 ratio difference vs ms1 intensity
plotXYDensity(log10(fData(sqaTMTSpectrum$eset)$ms1Int)
		, ratioDeltaSel
		, xlab="log10 MS1 intensity"
		, ylab="Compression Factor" )
abline(h=0, lty=2)


# plot log2 ratio difference vs ms1 interference
plotXYDensity(fData(sqaTMTSpectrum$eset)$interference
		, ratioDeltaSel
		, xlab="MS1 interference"
		, ylab="Compression Factor" )
abline(h=0, lty=2)

## plot minInt vs. interference
plotXYDensity(log(tmtInterference)[ tmtInterference != tmtMinInt],log(tmtMinInt)[ tmtInterference != tmtMinInt])
#abline(h=0, lty=2)


boxplot(ratioDeltaSel)
abline(h=0, lty=2)

plotVolcano(sqaTMTProtein, adjusted=T)
plotVolcano(sqaLFQProtein, adjusted=F)

print("DONE")

