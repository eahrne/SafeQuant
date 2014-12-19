source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/SafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Graphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")

library("affy")
library("limma")
library("gplots")

load("/Users/erikahrne/dev/R/workspace/SafeQuant/data/brozProeinCalMixTMT6.rda" )



#sp|P00722|BGAL_ECOLI	Beta-galactosidase OS=Escherichia coli (strain K12) GN=lacZ PE=1 SV=2	59	59	0	14.3898305085	116391.100660001	1024
#sp|P00489|PYGM_RABIT	Glycogen phosphorylase, muscle form OS=Oryctolagus cuniculus GN=PYGM PE=1 SV=3	65	65	51	10.6	97209.8471700004	843
#sp|P02789|TRFE_CHICK	Ovotransferrin OS=Gallus gallus PE=1 SV=2	50	50	3	12.1	77708.3879200003	705
#sp|P01012|OVAL_CHICK	Ovalbumin OS=Gallus gallus GN=SERPINB14 PE=1 SV=2	22	22	3	15.9545454545	42835.49196	386
#tr|B6V3I5|B6V3I5_BOVIN	Alpha-lactalbumin OS=Bos taurus PE=2 SV=1	7	7	0	12	16217.89797	142
#sp|P02666|CASB_BOVIN	Beta-casein OS=Bos taurus GN=CSN2 PE=1 SV=2	9	8	0	11.7777777778	25073.24305	224
#sp|P02662|CASA1_BOVIN	Alpha-S1-casein OS=Bos taurus GN=CSN1S1 PE=1 SV=2	13	12	0	11.9230769231	24495.41514	214
#sp|P02663|CASA2_BOVIN	

calMixProteinNames <- c(
		
		"sp|P00489|PYGM_RABIT"
		,"sp|P02789|TRFE_CHICK"
		,"sp|P01012|OVAL_CHICK"
		,"sp|P02666|CASB_BOVIN"
		,"sp|P00722|BGAL_ECOLI"
		,"tr|B6V3I5|B6V3I5_BOVIN"

)

expectedRatios <- list()
expectedRatios[calMixProteinNames] <- c(4,0.25,2,0.5,2,0.5)

#expectedRatios["sp|P00489|PYGM_RABIT"]

fData(esetTMT6Spectrum)$isNormAnchor <- !(fData(esetTMT6Spectrum)$proteinName %in% calMixProteinNames)

esetTMT6SpectrumNorm <- normalize(esetTMT6Spectrum, method="global")

pdf("/Users/erikahrne/tmp/tmp.pdf")

.qcPlots(esetTMT6Spectrum )
.qcPlots(esetTMT6SpectrumNorm)

par(mfrow=c(2,2))

for(targetProtein in calMixProteinNames){
	
	c1 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),1]
	c2 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),2]
	c3 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),3]
	c4 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),4]
	c5 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),5]
	c6 <- exprs(esetTMT6SpectrumNorm)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),6]
	
	r1 <- c1 / c2 
	r2 <- c3 / c4 
	r3 <- c5 / c6 
	
	main <- paste(targetProtein,"\nN=",length(c1),sep="")
	boxplot(log2(exprs(esetTMT6Spectrum)[(fData(esetTMT6Spectrum)$proteinName %in% targetProtein),]), main=main)
	boxplot(log2(r1),log2(r2),log2(r3), ylim=c(-2,2), main=main)
	abline(h=log2(expectedRatios[[targetProtein]]))
	
}

par(mfrow=c(1,1))

dev.off()

print("DONE")
