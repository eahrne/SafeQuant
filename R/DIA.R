#' get swath sizes ensuring equal number of precurosrs per swath
#' @param mzs precucrsor mz values
#' @param nbSwaths default 30
#' @param lowerOverlap default 2 lower bound overlap preceeding window
#' @return data.frame "binMean","lower","upper","delta"
#' @export
#' @note  No note
#' @details No details
#' @references No ref
#' @examples print("No examples")
getSwaths = function(mzs,nbSwaths=30,lowerOverlap = 2){
  
  swathCoord = levels(cut2(mzs,g=nbSwaths, digits=4)) %>% gsub("(\\[)|(\\)|( )|(\\]))","",.) %>% strsplit(",") %>% data.frame %>% t() %>% apply(.,2,as.numeric) %>% data.frame
  colnames(swathCoord) = c("lower","upper")
  swathCoord = cbind(binMean=rowMeans(swathCoord),swathCoord )
  swathCoord$delta = swathCoord$binMean - swathCoord$lower
  
  # overlap lower bound
  swathCoord$adjustedLower = swathCoord$lower -lowerOverlap
  swathCoord$adjustedBinMean = rowMeans(cbind(swathCoord$adjustedLower,swathCoord$upper))
  swathCoord$adjustedDelta = swathCoord$delta+(lowerOverlap/2)
  swathCoord$adjustedUpper = swathCoord$adjustedBinMean + swathCoord$adjustedDelta
  
  return(swathCoord)
}

