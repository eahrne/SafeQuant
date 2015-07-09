# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testExpressionAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testGraphics.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testIdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testSafeQuantAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testTMT.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testParser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/tests/testUserOptions.R")

cat("\n ---------------------- ALL TESTS PASSED ------------------ \n")


head(exprs(eset))


exprs(eset)[1,] <- c(1,1,2,2,3,3)

pData(eset)$isControl <- c(T,T,F,F,F,F) 

r <- getRatios(eset)

head(r)