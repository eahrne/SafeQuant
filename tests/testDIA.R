### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
  setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END


### TEST FUNCTIONS

tesGetSwaths <- function(){
  
  cat("--- tesGetSwaths: --- \n")
  

  set.seed(123)
  mzs= rnorm(100,500,100)
  
  nbSwaths = 4
  swaths = getSwaths(mzs,nbSwaths = nbSwaths, lowerOverlap=2 )
  
  hUnAdjusted = hist(mzs,breaks= c(min(mzs), swaths$lower[2:nbSwaths] , max(mzs)),freq = T )
  
  all(hUnAdjusted$counts >= 24) %>% stopifnot
  
  # plot
  # par(mfrow=c(3,1))
  # h = hist(mzs)
  # plot(hUnAdjusted, freq=T)
  # hist(mzs,breaks= c(swaths$adjustedLower , max(mzs)),freq = T )
  
  cat("--- tesGetSwaths: PASS ALL TEST --- \n")
  
  # count = c()
  # for(i in 1:nrow(swaths)){
  #   count = c(count,  sum( mzs >= swaths$lower[i] & mzs <= swaths$upper[i] ) )
  # }

}


### RUN TESTS

tesGetSwaths()