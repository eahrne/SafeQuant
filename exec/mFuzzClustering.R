# mFuzz cluster SafeQunat results
# 
# Author: erik ahrne
###############################################################################

### LOAD LIBRARIES
library(Mfuzz)

### LOAD LIBRARIES END

### LOAD DATA
sqRDataFile <- "/tmp//SQ_Results/SQ_Results_SQ.rData"
load(sqRDataFile)
eset <- sqaProtein$eset
# remove all objects we don't need
rm(list=ls()[ls() != "eset"])
### LOAD DATA END

### PARAMS
nbClusters <- 4
min.mem <- .5
pdfFile <- "/tmp/tmp.pdf"
xlsExportFolder <- "/tmp/"

# graphics
colo <- "fancy" 
ylab <- "Standardized Intensity"
ylim <- c(-1.5,1.5)
xlab <- "Condition" 

### PARAMS END

# standardise
# good idea when including all data?
eset <- standardise(eset)

### estimate fuzzification parameter using Schwammle and Jensen method
#       The algorithm needs a fuzzification parameter m in the range [1,n] which
#       determines the degree of fuzziness in the clusters. When
#       m reaches the value of 1 the algorithm works like a crisp
#       partitioning algorithm and for larger values of m the
#       overlapping of clusters is tend to be more.
m <- mestimate(eset)

### CLUSTER
#       This function is the core function for soft clustering. It groups genes based on the Euclidean distance
#       and the c-means objective function which is a weighted square error function. Each gene is assigned
#       a membership value between 0 and 1 for each cluster. Hence, genes can be assigned to different
#       clusters in a gradual manner. This contrasts hard clustering where each gene can belongs to a single
#       cluster.
#       Algorithm in brief:
#       1. Choose the number k of clusters.
#       2. Randomly generate k clusters and determine the cluster centers, or directly
#       generate k random points as cluster centers.
#       3. Assign each point to the nearest cluster center.
#       4. Recompute the new cluster centers.
#       5. Repeat the two previous steps until some convergence criterion is met (usually
#       that the assignment has not changed).
fuzzification <- mfuzz(eset,c=nbClusters,m=m)
### GRAPHICS
pdf(pdfFile)

### PLOT CLUSTERS WITH TREND LINE
par(mfrow=c(3,3))
for(clusterNb in 1:nbClusters){
	mfuzz.plot2(eset,cl=fuzzification, mfrow=NA, min.mem=min.mem, x11 = FALSE, ylim=ylim, single=clusterNb, ylab=ylab, xlab=xlab,colo=colo)
	lines(1:ncol(eset),fuzzification$centers[clusterNb,], lwd=2)
}
mfuzzColorBar(col=colo,main="Membership",cex.main=1)
dev.off()
cat("Figures exported to", pdfFile, "\n")
### GRAPHICS END

### XLS EXPORT
for(i in 1:nbClusters){
	
	xlsExportFile = paste(xlsExportFolder,"cluster",i,".xls" ,sep="")
	clusterMembers <- names(fuzzification$cluster[fuzzification$cluster == i])
	membershipValue <- fuzzification$membership[fuzzification$cluster == i,i]
	outDf <- data.frame(as.vector(unlist(membershipValue)),exprs(eset)[clusterMembers,])
	names(outDf)[1] <-  "membershipValue"
	outDf <- outDf[outDf[,1] >= min.mem,]
	write.table(file=xlsExportFile,outDf, row.names=F,sep="\t")
	cat("Created file", xlsExportFile, "\n")
}
### XLS EXPORT END
cat("DONE\n")