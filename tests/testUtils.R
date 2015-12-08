


### INIT
setwd(dirname(sys.frame(1)$ofile))
source("initTestSession.R")
### INIT END

testGetSQRootDir <- function(){
	
	cat(" --- testGetSQRootDir --- \n")
	stopifnot(dir.exists(getSQRootDir()))
	print(getSQRootDir())
	cat(" --- testGetSQRootDir: PASS ALL TEST --- \n")
	
}

testGetSQRootDir()


