# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### LOAD LIBRARIES

library(Mfuzz)
library(limma)
library(affy)

### LOAD LIBRARIES END

### SOURCE SCRIPTS

source('/Users/erikahrne/dev/R/workspace/SafeQuant/R/ExpressionAnalysis.R')

### SOURCE SCRIPTS END

### PARAMS

### filter
qvalueCutOff <- 0.05
log2RatioCutOff <- log2(1)

### mfuzz
nbClusters <- 8
min.mem <- .5
m <- 2

sqFile <- '/Users/erikahrne/dev/R/workspace/SafeQuant/inst/testData/mFuzzClustering/roland_20150205/data/sign_one_timepoint_qvalue0.01.tsv'
timeLabels = c(0,0.5,1,2,4,8)
pdfFile <- '/Users/erikahrne/tmp/mFuzz.pdf'
txtExportFolder <- '/Users/erikahrne/tmp/'
evalNbClusters <- T

### time course plots
colo <- "fancy" 
xlab <- "Time (min)" 
ylab <- "Norm Log Ratio"
ylim <- c(-1.5,1.5)

### PARAMS END

### PARSE DATA AND CREATE ESET

sqRes <- read.csv(sqFile,head=T,sep="\t")
peptideIdentifiers <-  paste(sqRes$peptide,sqRes$ptm,sep=":")
rownames(sqRes) <- peptideIdentifiers

qValues <- sqRes[,grepl("qValue",names(sqRes))]
rownames(qValues) <- peptideIdentifiers

ratios <- sqRes[,grepl("log2ratio",names(sqRes))]
rownames(ratios) <- peptideIdentifiers

### filter, keep only peptide d.e at at least one ime point
ok <- as.vector(unlist(apply((abs(ratios) > log2RatioCutOff ) & (qValues < qvalueCutOff),1,sum) > 0))

mfuzzInputMatrix <- as.matrix(data.frame(rep(0,nrow(ratios)) , ratios))
colnames(mfuzzInputMatrix) <- timeLabels
expDesign <- data.frame(isControl=c(T,rep(F,ncol(mfuzzInputMatrix)-1))
						,condition=paste("Cond",1:ncol(mfuzzInputMatrix))
						,row.names=colnames(mfuzzInputMatrix)
						)
				
eset <- createExpressionDataset(mfuzzInputMatrix[ok,],expDesign=expDesign,featureAnnotations=sqRes[ok,!grepl("log2ratio",names(sqRes))])	


### PARSE DATA AND CREATE ESET END


### MAIN

eset <- standardise(eset) ### STANDARDISE
### estimate using Schwammle and Jensen method
m <- mestimate(eset)
### cluster
fuzzification <- mfuzz(eset,c=nbClusters,m=m)




### MAIN

### MAIN END

### GRAPHICS

pdf(pdfFile)

### cluster size barplot
#barplot(fuzzification$size, names=1:nbClusters ,main="Cluster Sizes", ylab="# peptides")
barplot(apply((fuzzification$membership == apply(fuzzification$membership,1,max)) & (fuzzification$membership > min.mem),2,sum),names=1:nbClusters ,main="Cluster Sizes", ylab="# peptides")

#names(fuzzification)
#plot(fuzzification$centers)

### PLOT CLUSTERS WITH TREND LINE
par(mfrow=c(3,3))
for(clusterNb in 1:nbClusters){
	
	mfuzz.plot2(eset,cl=fuzzification, mfrow=NA, min.mem=min.mem, x11 = FALSE, ylim=ylim, single=clusterNb, ylab=ylab, xlab=xlab, timeLabels=timeLabels,colo=colo)
	#mfuzz.plot2(eset,cl=fuzzification, mfrow=NA, min.mem=min.mem, x11 = FALSE,single=clusterNb)
	lines(1:length(timeLabels),fuzzification$centers[clusterNb,], lwd=2)
}

mfuzzColorBar(col=colo,main="Membership",cex.main=1)

dev.off()
cat("Figures exported to", pdfFile, "\n")
### GRAPHICS END

### CSV EXPORT

### PER CLUSTER EXPORT

header <- names(sqRes)

for(i in 1:nbClusters){
	
	txtExportFile = paste(txtExportFolder,"cluster",i,".csv" ,sep="")
	clusterMembers <- names(fuzzification$cluster[fuzzification$cluster == i])
	membershipValue <- fuzzification$membership[fuzzification$cluster == i,i]
	
	outDf <- data.frame(as.vector(unlist(membershipValue)),sqRes[clusterMembers,])
	names(outDf)[1] <-  "membershipValue"
	
	outDf <- outDf[outDf[,1] >= min.mem,]
	
	write.table(file=txtExportFile,outDf, row.names=F,sep="\t")
}

cat("Results written to", txtExportFolder, "\n")
### CSV EXPORT END


####################### EVALUATE NB CLUSTERS
if(evalNbClusters){
	pdf(gsub("\\.pdf$","_CLUSTER_NB_EVAL.pdf",pdfFile))
	source("/Users/erikahrne/dev/R/workspace/SafeQuant/exec/mFuzzyClusteringClusterNbEval.R")
	dev.off()
}
####################### EVALUATE NB CLUSTERS END
cat("DONE\n")