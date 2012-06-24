# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### check if protein is contaminant
isCon <- function(accessions){
	return(regexpr("\\Wcon_",accessions,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("^con_",accessions,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("_con",accessions,ignore.case=TRUE, perl=TRUE) > -1 )
}

### check if protein is contaminant 
isDecoy <-function(accessions){
	
	return(((regexpr("rev_",accessions,ignore.case=TRUE) > -1)
						| (regexpr("decoy_",accessions,ignore.case=TRUE) > -1)
						| (regexpr("reverse_",accessions,ignore.case=TRUE) > -1)))
	
}

### calculate id level qvals based on target-decoy score distributions
getIdLevelQvals <- function(scores, isDecoy){
	
	targetScores = scores[!isDecoy]
	decoyScores = scores[isDecoy]
	
	qvals = c()
	for(score in scores){
		
		qval = sum(decoyScores >= score) / sum(targetScores >= score)
		qvals = c(qvals,qval)
	}
	
	return(qvals)
	
}
