# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
  setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END


### TEST FUNCTIONS

testCreateExpressionDataset <- function(){
  
  cat("--- testCreateExpressionDataset: --- \n")
  
  stopifnot(all.equal(pData(eset),expDesign))
  stopifnot(all.equal(sampleNames(eset),colnames(m)))
  stopifnot(all.equal(nrow(exprs(eset)),nrow(m)))
  cat("--- testCreateExpressionDataset: PASS ALL TEST --- \n")
  
}


testGetAllEBayes <- function(){
  
  
  cat("--- testGetAllEBayes: --- \n")
  p <- getAllEBayes(eset)
  stopifnot( mean(p[,"A"]) > mean(p[,"B"]) )
  p <- getAllEBayes(eset)
  stopifnot( mean(p[,"A"]) > mean(p[,"B"]) )
  pAdj <- getAllEBayes(eset, adjust=T)
  
  stopifnot( mean(pAdj[,"A"]) > mean(p[,"A"]) )
  
  ### paired desgn
  esetNonPaired <- esetPaired
  pData(esetNonPaired) <- pData(esetPaired)[,1:2]
  
  pValPaired <- getAllEBayes(esetPaired,adjust=F,method=c("all"))
  pValPairedPairwise <- getAllEBayes(esetPaired,adjust=F,method=c("pairwise"))
  
  pValNonPaired <- getAllEBayes(esetNonPaired,adjust=F,method=c("all"))
  pValNonPairedPairwise <- getAllEBayes(esetNonPaired,adjust=F,method=c("pairwise"))
  stopifnot(all(apply(pValPaired,2,sum) < apply(pValNonPaired,2,sum)))
  stopifnot(all(apply(pValPairedPairwise,2,sum) < apply(pValNonPairedPairwise,2,sum)))
  
  # adjustment filter
  
  adjustfilter <- data.frame(rep(T,nrow(eset)),rep(T,nrow(eset)) )
  adjustfilter[1,1] <- F
  adjustfilter[2,2] <- F
  
  pTmp <- getAllEBayes(eset,adjust=T, adjustFilter=adjustfilter)
  stopifnot(p[1,1] == pTmp[1,1])
  stopifnot(p[2,2] == pTmp[2,2])
  stopifnot(is.na(pTmp[3,2]))
  cat("--- testGetAllEBayes: PASS ALL TEST --- \n")
  
  # Question: limma multiple groups comparison produces different pvalue comparing with two group comparsion	
  # https://support.bioconductor.org/p/44216/	
  # https://support.bioconductor.org/p/60556/
}

testGetRatios <- function(){
  
  cat("--- testGetRatios: --- \n")
  
  r <- getRatios(eset)
  
  exprs(eset)[1,c(5:6,1:2)]
  #  C_rep_1  C_rep_2  A_rep_1  A_rep_2 
  # 9.963852 9.966003 9.965568 9.967214 
  
  ## NOTE: C is control
  stopifnot(all.equal(r[1,1],median(log2(exprs(eset))[1,1:2]) - median(log2(exprs(eset))[1,5:6]) ) )
  
  stopifnot(mean(r[,"A"]) < mean(r[,"B"]))
  stopifnot(all.equal(r,getRatios(eset,method="mean")))
  r <- getRatios(eset)
  stopifnot(all.equal(mean(r[,"B"]),mean(r[,"B"])))
  
  stopifnot(all(round(apply(getRatios(esetPaired, log2=F),2,mean),1) == c(1.5,1.2)))
  stopifnot(ncol(r)  == 2)
  
  rPaired <- getRatios(esetPaired, method="paired")
  stopifnot(ncol(rPaired)  == 4)
  stopifnot(all(	log2(exprs(esetPaired)[,1]) - log2(exprs(esetPaired)[,5]) == rPaired[,1] ))
  stopifnot(all(	log2(exprs(esetPaired)[,4]) - log2(exprs(esetPaired)[,6]) == rPaired[,4] ))
  stopifnot(all(	exprs(esetPaired)[,4] / exprs(esetPaired)[,6] == getRatios(esetPaired, method="paired", log2=F)[,4] ))
  # should fail	
  stopifnot((inherits(try(	getRatios(eset, method="paired"), silent = TRUE), "try-error")))
  
  
  stopifnot(all(round(getRatios(esetPaired)[,1],2) ==  round(apply(rPaired[,1:2],1,median),2)))
  
  
  cat("--- testGetRatios: PASS ALL TEST --- \n")
  
}

testGetAllCV <- function(){
  
  cat("--- testGetAllCV: --- \n")
  
  cv <- getAllCV(eset)
  
  stopifnot(all.equal(cv[1,"A"] , (sd(exprs(eset)[1,unlist(pData(eset)$condition == "A")]) / mean(exprs(eset)[1,unlist(pData(eset)$condition == "A")]))))
  stopifnot(all.equal(cv[200,"C"] , sd(exprs(eset)[200,unlist(pData(eset)$condition == "C")]) / mean(exprs(eset)[200,unlist(pData(eset)$condition == "C")])))
  
  cvPaired <-  getAllCV(esetPaired)
  stopifnot(all(is.na(cvPaired[,1])))
  stopifnot(all(!is.na(cvPaired[,2])))
  stopifnot(cvPaired[3,"A"] ==  (sd(exprs(esetPaired)[3,1:2] / exprs(esetPaired)[3,5:6]) /	mean(exprs(esetPaired)[3,1:2] / exprs(esetPaired)[3,5:6])))
  
  cat("--- testGetAllCV: PASS ALL TEST --- \n")
}

testGlobalNormalize <- function(){
  
  cat("--- testGlobalNormalize: --- \n")
  
  globalNormFactors <- getGlobalNormFactors(eset,method="sum")
  ### add normalization factors to ExpressionSet
  pData(eset) <- cbind(pData(eset),globalNormFactors)
  esetNorm <- globalNormalize(eset,globalNormFactors)
  
  stopifnot(all.equal(as.vector(unlist(apply(exprs(esetNorm),2,sum))),as.vector(rev(unlist(apply(exprs(esetNorm),2,sum))))))
  
  stopifnot(pData(esetNorm)$normFactor[1] == 1)
  stopifnot(pData(esetNorm)$normFactor[2] != 1)
  
  cat("--- testGlobalNormalize: PASS ALL TEST --- \n")
  
  # Should generate error and stop	
  #	fData(eset)$isNormAnchor <- rep(F,nrow(eset))
  #	esetNorm <- normalizeIntensities(eset)
  
}

testGetSignalPerCondition <- function(){
  
  cat("--- testGetSignalPerCondition: --- \n")
  stopifnot(sum(getSignalPerCondition(eset,method="min")[,"A"] <= getSignalPerCondition(eset,method="max")[,"A"]) == nrow(eset))
  stopifnot(sum(getSignalPerCondition(eset,method="min")[,"C"] <= getSignalPerCondition(eset,method="median")[,"C"]) == nrow(eset))
  stopifnot(getSignalPerCondition(eset,method="sd")[1,2] == sd(exprs(eset)[1,3:4]))
  cat("--- testGetSignalPerCondition: PASS ALL TEST --- \n")
}

testBaselineIntensity <- function(){
  
  cat("--- testBaselineIntensity: --- \n")
  
  allInt <- as.vector(unlist(exprs(eset)))
  bl <- round(getBaselineIntensity(allInt,promille=5),2)
  #hist(allInt)
  #abline(v=bl,lwd=2)
  stopifnot(all.equal(bl[[1]] , 997.82) )
  
  
  set.seed(1234)
  suppressWarnings(allInt2 <- log10(rnorm(1000,1,1)))
  bl2 <- round(getBaselineIntensity(allInt2,promille=5),2)
  stopifnot(all.equal(-1.95,bl2[[1]]))
  #hist(allInt2)
  #abline(v=bl2,lwd=2)
  
  cat("--- testBaselineIntensity: PASS ALL TEST --- \n")
  
}

testRollUp <- function(){
  
  cat(" --- testRollUp --- \n")
  
  uProt = unique( fData(eset)$proteinName )
  rollUpEset1 <- rollUp(eset,featureDataColumnName= c("proteinName"), method=c("sum"))
  stopifnot( length( uProt ) == nrow(rollUpEset1)) 
  t1 = subset(exprs(eset),fData(eset)$proteinName %in% uProt[115] ) %>% colSums
  stopifnot( all(t1 == subset(exprs(rollUpEset1),fData(rollUpEset1)$proteinName %in% uProt[115] )))
  t2 = subset(exprs(eset),fData(eset)$proteinName %in% uProt[23] ) %>% colSums
  stopifnot( all(t2 == subset(exprs(rollUpEset1),fData(rollUpEset1)$proteinName %in% uProt[23] )))
  
  rollUpEset2 <- rollUp(eset ,featureDataColumnName= c("ptm"), method=c("sum"))
  stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset2)) 
  
  rollUpEset3 <- rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("ptm"), method=c("mean"))
  stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset3)) 

  #print(exprs(rollUpEset2))
  
  stopifnot(all.equal(sum(exprs(rollUpEset2)),sum(exprs(rollUpEset1)))) ### test sum
  stopifnot(sum(exprs(rollUpEset1)) != sum(exprs(rollUpEset3))) ### test mean
  
  rollUpEset4 <- rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top3"))
  stopifnot(sum(exprs(rollUpEset3)) != sum(exprs(rollUpEset4))  ) ### test top 3
  
  esetTmp <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementCsvFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementCsvFile1))
  uProt3 = unique( fData(esetTmp)$proteinName )
  t3 = subset(exprs(esetTmp),fData(esetTmp)$proteinName %in% uProt3[2] ) %>% colSums
  stopifnot( all(t3 == subset(exprs(esetTmp),fData(esetTmp)$proteinName %in% uProt[2] )))
  
  
  # test rollup of NA_IMP
  stopifnot(ncol(getImputedIntensities(rollUp(rollUp(eset ,featureDataColumnName= c("peptide","ptm"), method=c("sum"))))) == ncol(eset))
  esetImp = eset
  exprs(esetImp)[1:20,1] = NA
  exprs(esetImp)[1:20,3] = NA
  exprs(esetImp)[21:30,2] = NA
  esetGMin = sqImpute(esetImp, method="gmin")
  rolledEsetGMin = rollUp(esetGMin)
  
  impInt = getImputedIntensities(esetGMin)
  impIntRolled = getImputedIntensities(rolledEsetGMin)
  stopifnot(all(colSums(impInt) == colSums(impIntRolled))) 

    # check single proteins
  stopifnot(impIntRolled[match(uProt[6], rownames(impIntRolled) ),] ==   colSums(impInt[fData(esetGMin)$proteinName %in% uProt[6]  ,]))
  stopifnot(impIntRolled[match(uProt[16], rownames(impIntRolled) ),] ==   colSums(impInt[fData(esetGMin)$proteinName %in% uProt[16]  ,]))
  
  cat(" --- testRollUp: PASS ALL TEST  --- \n")
  
}

testTopX <- function(){
  
  cat(" --- testTopX --- \n")
  
  entryData1  <- data.frame(t(matrix(c(1,1,1,3,3,3,2,2,2,5,5,5),ncol=4)))
  rownames(entryData1) <- paste("peptide",1:nrow(entryData1),sep="_")
  #           X1 X2 X3
  # peptide_1  1  1  1
  # peptide_2  3  3  3
  # peptide_3  2  2  2
  # peptide_4  5  5  5
  stopifnot(all.equal(rep(10/3,3) , as.vector(unlist(getTopX(entryData1))) ))
  
  entryData2 <- data.frame(t(matrix(c(1,1,1,3,3,3,2,2,2,5,5,NA),ncol=4)))
  rownames(entryData2) <- paste("peptide",1:nrow(entryData2),sep="_")
  
  #           X1 X2 X3
  # peptide_1  1  1  1
  # peptide_2  3  3  3
  # peptide_3  2  2  2
  # peptide_4  5  5 NA
  stopifnot(all.equal(c(4,4,3) ,  as.vector(unlist(getTopX(entryData2,topX=2)))))
  
  # 1 row
  stopifnot(all.equal(rep(1,3),as.vector(unlist(getTopX(entryData1[1,])))))
  
  # 1 col
  stopifnot(all.equal(getTopX(entryData1)[[1]],getTopX(entryData1[,1])))
  
  top1 <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top1"))),1,sum)
  top3 <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top3"))),1,sum)
  meanInt <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("mean"))),1,sum)
  
  stopifnot(sum(top1 >=  top3 ) == length(top3))
  stopifnot(sum(top1) > sum(top3))
  stopifnot(sum(top1) > sum(meanInt))
  stopifnot(all.equal(sum(top3),sum(meanInt)))
  
  cat(" --- testTopX: PASS ALL TEST  --- \n")
  
}

testGetIBAQEset <- function(){
  
  cat(" --- testGetIBAQEset --- \n")
  
  # read protein fasta
  proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
  
  iBaqEset <- getIBAQEset(eset,proteinDB=proteinDB)
  stopifnot(round(exprs(iBaqEset))[1,1] == 125)
  stopifnot(round(exprs(iBaqEset))[2,2] == 38)
  
  cat(" --- testGetIBAQEset: PASS ALL TEST --- \n")
}


testGetLoocvFoldError <- function(){
  
  cat(" --- testGetLoocvFoldError --- \n")
  #plotCalibrationCurve(fit)
  stopifnot(  round(sum(getLoocvFoldError(absEstSimData))) == -8)
  cat(" --- testGetLoocvFoldError: PASS ALL TEST --- \n")
}

testSqNormalize <- function(){
  
  cat(" --- testSqNormalize --- \n")
  
  stopifnot(nrow(sqNormalize(eset, method = "global")) == 900)
  esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
  stopifnot(nrow(sqNormalize(esetTmp, method = "rt")) == 496)
  
  cat(" --- testSqNormalize: PASS ALL TEST --- \n")
}


testRtNormalize <- function(){
  
  
  cat(" --- testRTNormalize --- \n")
  esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
  rtNormFactors <- getRTNormFactors(esetTmp, minFeaturesPerBin=100)
  stopifnot(nrow(rtNormalize(esetTmp,rtNormFactors)) == 97)
  
  # Stop if rtNormFactors doesn't cover all retention times.
  #rtNormalize(esetTmp,rtNormFactors[1:2,])
  
  cat(" --- testRTNormalize: PASS ALL TEST --- \n")
}

testRemoveOutliers <- function(){
  
  cat(" --- testRemoveOutliers --- \n")
  set.seed(1234)
  stopifnot(sum(is.na(removeOutliers(c(-10,  rnorm(100), 10)))) == 2)
  cat(" --- testRemoveOutliers: PASS ALL TEST --- \n")
}

testPerFeatureNormalization <- function(){
  
  cat(" --- testPerFeatureNormalization: --- \n")
  
  normFactors <- exprs(eset)[1:10,1:3]
  colnames(normFactors) <- c("A","B","C")
  rownames(normFactors) <- fData(eset)[1:10,]$proteinName
  normFactors[is.finite(normFactors)] <- 1
  eNorm <-  perFeatureNormalization(eset,normFactors)
  
  stopifnot(all.equal(exprs(eNorm)[1,1] , (exprs(eset)[1,1]-1)))
  stopifnot(all.equal(exprs(eNorm)[1,2] , (exprs(eset)[1,2]-1)))
  stopifnot(all.equal(exprs(eNorm)[13,2] , (exprs(eset)[13,2])))
  
  normFactors[,2] <- 200
  eNorm <-  perFeatureNormalization(eset,normFactors)
  stopifnot(all.equal(exprs(eNorm)[1,1] , (exprs(eset)[1,1]-1)))
  stopifnot(all.equal(exprs(eNorm)[1,3] ,  (exprs(eset)[1,3]-200)))
  stopifnot(all.equal(exprs(eNorm)[1,4] ,  (exprs(eset)[1,4]-200)))
  stopifnot(all.equal(exprs(eNorm)[1,5] ,  (exprs(eset)[1,5]-1)))
  
  normFactors[,3] <- 0
  eNorm <-  perFeatureNormalization(eset,normFactors)
  stopifnot(all.equal(exprs(eNorm)[1,5] , exprs(eset)[1,5]))
  stopifnot(all.equal(exprs(eNorm)[2,6] , exprs(eset)[2,6]))
  
  normFactors[3,] <- 1000 
  eNorm <-  perFeatureNormalization(eset,normFactors)
  stopifnot(all.equal(exprs(eNorm)[3,] , (exprs(eset)[3,]-1000)))
  stopifnot(all.equal(exprs(eNorm)[20:50,] , exprs(eset)[20:50,]))
  
  # re-order normFactors columns
  eNorm <-  perFeatureNormalization(eset,normFactors[,rev(colnames(normFactors))])
  stopifnot(all.equal(exprs(eNorm)[3,] , (exprs(eset)[3,]-1000)))
  stopifnot(all.equal(exprs(eNorm)[20:50,] , exprs(eset)[20:50,]))
  
  cat(" --- testPerFeatureNormalization: PASS ALL TEST --- \n")
  
  #coveredPeptideSel <- fData(eset)$proteinName %in% rownames(normFactors)
  #exprs(eset)[coveredPeptideSel,]	<- exprs(eset)[coveredPeptideSel, ] - normFactors[as.character(fData(eset)[coveredPeptideSel,]$proteinName),pData(eset)$condition]
  
}


testGetMaxIndex <-function(){
  
  cat(" --- testGetMaxIndex:  --- \n")
  d <- data.frame(s=c(NA,NA,NA,NA,1,1:4),lab=sort(rep(c("A","B","C"),3)))
  DT <- data.table(d)
  setkey(DT,lab)
  
  
  stopifnot(all.equal(c(1,5,9) , DT[, .I[getMaxIndex(s)], by=lab ]$V1  ))
  cat(" --- testGetMaxIndex: PASS ALL TEST  --- \n")
  
}

testCreatePairedExpDesign <- function(){
  
  cat(" --- testCreatePairedExpDesign:  --- \n")
  stopifnot(all( pData(createPairedExpDesign(eset))$subject == as.factor(rep(1:2,3))))
  stopifnot((inherits(try(createPairedExpDesign(eset[,1:3]), silent = TRUE), "try-error")))
  stopifnot((inherits(try(createPairedExpDesign(eset[,1:4]), silent = TRUE), "ExpressionSet")))
  cat(" --- testCreatePairedExpDesign: PASS ALL TEST  --- \n")
  
}

testGetFTestPValue <- function(){
  
  cat(" --- testGetFTestPValue:  --- \n")
  
  set.seed(1234)
  # creat test dataset
  nbFeatures <- 100
  AA1 <-rnorm(nbFeatures,100,100*0.1)
  AA2 <-rnorm(nbFeatures,100,100*0.1)
  BB1 <- rnorm(nbFeatures,100,100*0.1)
  BB2 <- rnorm(nbFeatures,100,100*0.1)
  CC1 <-rnorm(nbFeatures,100,100*0.1)
  CC2 <- rnorm(nbFeatures,100,100*0.1)
  AA1[1:10] <- AA1[1:10]*2
  AA2[1:10] <- AA2[1:10]*2
  BB1[11:20] <- BB1[11:20]*3
  BB2[11:20] <- BB2[11:20]*3
  
  mF <- log10(as.matrix(data.frame(AA1,AA2,BB1,BB2,CC1,CC2) ))
  colnames(mF) <- c("A_rep_1","A_rep_2","B_rep_1","B_rep_2","C_rep_1","C_rep_2")
  expDesign <- data.frame(condition=c("A","A","B","B","C","C"),isControl=c(T,T,F,F,F,F),row.names=colnames(m))
  
  featureAnnotations <- data.frame(row.names=as.character(1:nbFeatures), bla=as.character(1:nbFeatures), blabla=1:nbFeatures)
  esetF <- createExpressionDataset(expressionMatrix=mF,expDesign=expDesign,featureAnnotations=featureAnnotations)
  
  pValueF = getFTestPValue(esetF, log=F, adjust=F)
  
  stopifnot(!all(pValueF[21:100] < 0.001))
  stopifnot(all(pValueF[1:20] < 0.001))
  
  # plot dataset
  #gplots::heatmap.2(exprs(esetF), dendrogram = "row",trace = "none")
  
  # low p-values should be assigned 1:20
  #plot( pValueF , log="y", col=c( rep("red",10),rep("blue",10),rep("grey",80)), pch =19, ylab ="p-value" )
  
  cat(" --- testGetFTestPValue: PASS ALL TEST  --- \n")
  
}

testRollUpDT <- function(){
  
  cat(" --- testRollUpDT --- \n")
  
  rollUpEset1 <- rollUpDT(eset,featureDataColumnName= c("proteinName"), method=c("sum"))
  stopifnot( length( unique( fData(eset)$proteinName ) ) == nrow(rollUpEset1)) 
  
  rollUpEset2 <- rollUpDT(eset ,featureDataColumnName= c("ptm"), method=c("sum"))
  stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset2)) 
  
  rollUpEset3 <- rollUpDT(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("ptm"), method=c("mean"))
  stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset3)) 
  
  #print(exprs(rollUpEset2))
  
  stopifnot(all.equal(sum(exprs(rollUpEset2)),sum(exprs(rollUpEset1)))) ### test sum
  stopifnot(sum(exprs(rollUpEset1)) != sum(exprs(rollUpEset3))) ### test mean
  
  rollUpEset4 <- rollUpDT(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top3"))
  stopifnot(sum(exprs(rollUpEset3)) != sum(exprs(rollUpEset4))  ) ### test top 3
  
  cat(" --- testRollUpDT: PASS ALL TEST  --- \n")
  
  # esetTmp <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementCsvFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementCsvFile1))
  # rollUpEsetProteinAllAccessions <- rollUpDT(esetTmp,featureDataColumnName= c("proteinName"), method=c("sum"))
  # stopifnot(fData(rollUpEsetProteinAllAccessions)$allAccessionsTMP[68] == "sp|E9PAV3|NACAM_HUMAN;sp|Q9BZK3|NACP1_HUMAN")
  # 
  # 
  # rollUpEsetProteinAllAccessions <- rollUpDT(esetTmp,featureDataColumnName= c("proteinName"), method=c("sum"))
  # rollUpDTEsetProteinAllAccessions <- rollUpDT(esetTmp,featureDataColumnName= c("proteinName"), method=c("sum"))
  
  
}

testSqImpute <-function(){
  
  esetImp = eset

  exprs(esetImp)[1:20,1] = NA
  exprs(esetImp)[1:20,3] = NA
  exprs(esetImp)[21:30,2] = NA

  cat(" --- testSqImpute:  --- \n")
  esetGMin = sqImpute(esetImp, method="gmin")
  esetLMin = sqImpute(esetImp, method="lmin")
  esetKNN = sqImpute(esetImp, method="knn")
  esetPPCA = sqImpute(esetImp, method="ppca")
  esetGMean = sqImpute(esetImp, method="gmean")
  esetLMean = sqImpute(esetImp, method="lmean")
  
  stopifnot(sum(is.na(exprs(esetGMin))) == 0)
  stopifnot(sum(is.na(exprs(esetLMin))) == 0)
  stopifnot(sum(is.na(exprs(esetKNN))) == 0)
  stopifnot(sum(is.na(exprs(esetPPCA))) == 0)
  stopifnot(sum(is.na(exprs(esetGMean))) == 0)
  stopifnot(sum(is.na(exprs(esetLMean))) == 0)
  
  # fData(esetGMin)[ grepl("^NA\\_IMP",names(fData(esetGMin)))  ] %>% head
  # fData(esetGMean)[ grepl("^NA\\_IMP",names(fData(esetGMean)))  ] %>% head
  # fData(esetKNN)[ grepl("^NA\\_IMP",names(fData(esetKNN)))  ] %>% head
  # fData(esetPPCA)[ grepl("^NA\\_IMP",names(fData(esetPPCA)))  ]  %>% head
  #
  

  
  cat(" --- testSqImpute: PASS ALL TEST  --- \n")
  
}




testGetImputedIntensities = function(){
  
  cat(" --- testGetImputedIntensities: --- \n")
  
  # no imputed values
  stopifnot(all(getImputedIntensities(eset) == 0 ))
  
  # with imputed values
  
  esetImp = eset
  exprs(esetImp)[1:20,1] = NA
  exprs(esetImp)[1:20,3] = NA
  exprs(esetImp)[21:30,2] = NA
  esetGMin = sqImpute(esetImp, method="gmin")
  isNA = is.na(exprs(esetImp))
  mImp =   getImputedIntensities(esetGMin)
  stopifnot(all(mImp[isNA] > 0))
  stopifnot(all(mImp[!isNA] == 0))
  
  cat(" --- testGetImputedIntensities: PASS ALL TEST  --- \n")
  
}

testGetNAFraction = function(){
  
  cat(" --- testGetNAFraction:  --- \n")
  
  esetImp = eset
  exprs(esetImp)[1:20,1] = NA
  exprs(esetImp)[1:20,3] = NA
  exprs(esetImp)[21:30,2] = NA
  esetGMin = sqImpute(esetImp, method="gmin")
  
  fracRun = getNAFraction(esetGMin, method="run") 
  fracCond =  getNAFraction(esetGMin, method="cond")
  fracRatio =  getNAFraction(esetGMin) 
  
  stopifnot(ncol(fracRun) == 6 )
  stopifnot(ncol(fracCond) == 3 )
  stopifnot(ncol(fracRatio) == 2 )
  
  # run
  stopifnot(all(colSums(fracRun[1:30,1:3]) == c(20,10,20)))
  # ratio
  stopifnot(sum(exprs(esetGMin)[3,c(1)]) / sum(exprs(esetGMin)[3,c(1,2,5,6)]) == fracRatio[3,1])
  # cond
  stopifnot(sum(exprs(esetGMin)[20,c(3)]) / sum(exprs(esetGMin)[20,3:4]) == fracCond[20,2]) 
  
  cat(" --- testGetNAFraction: PASS ALL TEST  --- \n")
  
}

testGetNAFraction()

### TEST FUNCTIONS END

#### TESTS
#testGetRatios()
#testGetAllEBayes()
#testGetSignalPerCondition()
#testGetRatios()

#testRollUp()


if(F){
  testCreateExpressionDataset()
  testGetAllEBayes()
  testGetRatios()
  testGetAllCV()
  testGlobalNormalize()
  #testNormalise()
  testRtNormalize()
  
  testBaselineIntensity()
  testRollUp()
  testRollUpDT()
  
  testTopX()
  
  testGetLoocvFoldError()
  testRemoveOutliers()
  testPerFeatureNormalization()
  testGetMaxIndex()
  testGetIBAQEset()
  testCreatePairedExpDesign()
  testGetFTestPValue()
  
  testSqImpute()
  testGetImputedIntensities()
  
}

### get fraction missing values per peptide and protein.





