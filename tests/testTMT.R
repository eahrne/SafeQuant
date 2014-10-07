# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/Parser.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/TMT.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/IdentificationAnalysis.R")
source("/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R")

### test functions

testGetImpuritiesMatrix <- function(){
	cat(" --- testGetImpuritiesMatrix --- \n")
	# 6-plex
	# old stopifnot(0.094 ==  getImpuritiesMatrix(6)[1,2])
	stopifnot(0.06 ==  round(getImpuritiesMatrix(6)[1,2],2))
	# 10-plex
	stopifnot(0.047 ==  round(getImpuritiesMatrix(10)[1,2],3))
	cat(" --- testGetImpuritiesMatrix: PASS ALL TEST --- \n")
}

testPurityCorrectTMT <- function(){

	cat(" --- testPurityCorrectTMT --- \n")
	# 6-plex
	# old stopifnot( round(9.998839,4)  ==  round(purityCorrectTMT(tmtTestData6Plex,impurityMatrix=getImpuritiesMatrix(6))[1,1],4))
	stopifnot( round(8.999833,4)  ==  round(purityCorrectTMT(tmtTestData6Plex,impurityMatrix=getImpuritiesMatrix(6))[2,1],4))
	
	# 10-plex
	stopifnot( 10  ==  round(purityCorrectTMT(tmtTestData10Plex,impurityMatrix=getImpuritiesMatrix(10))[1,1],4))
	cat(" --- testPurityCorrectTMT: PASS ALL TEST --- \n")
	
}

#testGetIntSumPerProtein <- function(){
#	
#	cat(" --- testGetIntSumPerProtein --- \n")
#	SCAFFOLDOUTDATA <- parseScaffoldRawFile(tmt6PlexRawTestFile)
#	UNGROUPEDRAWINTDATA <- SCAFFOLDOUTDATA[,10:15]
#	proteinSummary <- getIntSumPerProtein(UNGROUPEDRAWINTDATA,SCAFFOLDOUTDATA$Accession.Numbers,SCAFFOLDOUTDATA$Peptide.Sequence)
#	stopifnot(proteinSummary$peptidesPerProtein[3] == 43) 
#	cat(" --- testGetIntSumPerProtein: PASS ALL TEST --- \n")	
#}

testCreateExpDesign <- function(){
	
	cat(" --- testCreateExpDesign --- \n")

	stopifnot(sum(createExpDesign("1,2,3:4,5,6",6)$isControl == c(T,T,T,F,F,F)) == 6 )
	stopifnot(sum(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$isControl == c(T,T,F,F,F,F,F,F,F,F)) == 10 )
	stopifnot(sum(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$condition[1:2] == c("Ctrl","Ctrl")) == 2 )
	stopifnot(length(unique(createExpDesign("10,2:3:4,5:6,7,8:9,1",10)$condition)) == 5 )
	
	#createExpDesign("1,4,10:2,5,8:3,6,9",10)
	cat(" --- testCreateExpDesign: PASS ALL TEST --- \n")
	
}

### compare our impurity correction to MSnbase
comparePurityCorrectionToMsnbase <- function(){
	
	cat(" --- comparePurityCorrectionToMsnbase --- \n")
	
	library(MSnbase)
	impurities <- matrix(c(0.929, 0.059, 0.002, 0.000,
					0.020, 0.923, 0.056, 0.001,
					0.000, 0.030, 0.924, 0.045,
					0.000, 0.001, 0.040, 0.923),
			nrow = 4)
	
	qnt <- quantify(itraqdata[1:2,],reporters = iTRAQ4)
	
	exprs(qnt)[1,] <- rep(1,4)
	exprs(qnt)[2,] <- 1:4
	qnt.crct <- purityCorrect(qnt, impurities)
	
	### our correction method
	testA <- as.vector(t(solve(impurities) %*% rep(1,4)))
	testB <- as.vector(t(solve(impurities) %*% 1:4))
	
	stopifnot( sum(round(testA,3) == round(exprs(qnt.crct)[1,],3)) == 4 )
	stopifnot( sum(round(testB,3) == round(exprs(qnt.crct)[2,],3)) == 4 )
	
	cat(" --- comparePurityCorrectionToMsnbase: PASS ALL TEST --- \n")
	
} 

### test functions end

# INIT

### CREATE TEST DATA

tmtTestData6Plex <- matrix(rep(10,24),ncol=6)
tmtTestData6Plex[2,1:3] <- c(9,9,9) 
tmtTestData6Plex[3,1:3] <- c(100,100,100) 
tmtTestData6Plex[4,c(1,3,5)] <- c(100,100,100) 

tmtTestData10Plex <- matrix(rep(10,100),ncol=10)

### CREATE TEST DATA END

testDir <- dirname(sys.frame(1)$ofile)
testDir <- gsub("tests\\/tmp","inst/tests/",testDir)
tmt6PlexRawTestFile <- paste(testDir,"/TMT_6-Plex_Scaffold_Raw_Export_Example.xls",sep="")
tmt10PlexRawTestFile <- paste(testDir,"/TMT_10-Plex_Scaffold_Raw_Export_Example.xls",sep="")

### INIT END

### TESTS

testGetImpuritiesMatrix()
testPurityCorrectTMT()
#testGetIntSumPerProtein()
testCreateExpDesign()

#comparePurityCorrectionToMsnbase()
### TESTS END






