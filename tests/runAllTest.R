# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### Don't run if in CHECK mode
if(!grepl("SafeQuant\\.Rcheck$",getwd())){
	
	setwd(dirname(sys.frame(1)$ofile))
	source("initTestSession.R")
	source("testExpressionAnalysis.R")
	source("testGraphics.R")
	source("testIdentificationAnalysis.R")
	source("testSafeQuantAnalysis.R")
	source("testTMT.R")
	source("testParser.R")
	source("testUserOptions.R")
	
	cat("\n ---------------------- ALL TESTS PASSED ------------------ \n")
	
}






