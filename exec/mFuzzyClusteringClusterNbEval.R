# TODO: Add comment
# 
# Author: erikahrne
###############################################################################


####################### EVALUATE NB CLUSTERS

nbClusterRange <- 3:12
rounds <- 3

### within cluster error
wErrorPerClusterNb <- data.frame(row.names=1:rounds)
for(i in nbClusterRange){
	wError <- c()
	for(j in 1:rounds){
		wErrorLoc <- mfuzz(eset,c=i,m=m)$withinerror
		wError <- c(wError,wErrorLoc )
	}
	wErrorPerClusterNb <- data.frame(wErrorPerClusterNb,wError)
}


par(mfrow=c(1,3))

### plot within cluster error
names(wErrorPerClusterNb) <- nbClusterRange
boxplot(wErrorPerClusterNb, main="min J(c,m), within cluster error", ylab="within cluster error")

### empty clusters for a given number of selected clusterss
cselection(eset,m,crange=nbClusterRange,repeats=5,visu=TRUE)

#The average minimum centroid distance for the given range of cluster
Dmin(eset,m,crange=nbClusterRange,repeats=rounds,visu=TRUE)

if(T){
	### various evaluation indices
#Calculates the values of several fuzzy validity measures. The values of the indexes can be independently used in order to evaluate and compare clustering partitions or even to determine the number
#of clusters existing in a data set.

	idxName <- names(fclustIndex(mfuzz(eset[1:100,],c=3,m=m),mfuzzInputMatrix[1:100,]))
	fclustIndexOut <- data.frame(row.names=idxName)
		
	for(nbCl in nbClusterRange){
	
		calculatedIndices <- fclustIndex(mfuzz(eset,c=nbCl,m=m),mfuzzInputMatrix)
		fclustIndexOut <- data.frame(fclustIndexOut,calculatedIndices)
		
	}
	
	fclustIndexOut <- data.frame(t(as.matrix(fclustIndexOut)))
	#names(fclustIndexOut) <- idxNames
	
	par(mfrow=c(2,4))
	
	for(i in 1:ncol(fclustIndexOut)){
		
	 	lab <- idxName[i]
		
	    index <- fclustIndexOut[,i]
		
		ok <- is.finite(index)
		
		if(sum(ok) > 0){
			plot(nbClusterRange,fclustIndexOut[,i], main=lab, ylab=lab)
		}else{
			plot(nbClusterRange,rep(0,length(nbClusterRange)), main=lab, ylab=lab)
		}
	}
}

### various evaluation indices END

####################### EVALUATE NB CLUSTERS END
