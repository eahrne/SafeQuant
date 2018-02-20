# TODO: Add comment
# 
# Author: eahrne
###############################################################################



### Don't run if in CHECK mode
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	
	setwd(dirname(sys.frame(1)$ofile))
	source("initTestSession.R")
	source("testParser.R")
	source("testIdentificationAnalysis.R")
	source("testExpressionAnalysis.R")
	source("testSafeQuantAnalysis.R")
  source("testTMT.R")
  source("testUserOptions.R")
	source("testGraphics.R")
	source("testTargeted.R")
  source("testGGGraphics.R")
  source("testDIA.R")
	
	cat("\n ---------------------- ALL TESTS PASSED ------------------ \n")
	
}






