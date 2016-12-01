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

#install.packages("/Users/erikahrne/dev/R/workspace/SafeQuant/", repos = NULL, type="source")
#library(SafeQuant)


### INIT
skylineExportFile <- "testData/skyline_transition_results.csv"
df = read.csv(skylineExportFile)
df$sample = substr(df$Replicate.Name,1,3)

### create transition  eset

# list all transitions "y10:1:IADFGLATQLK"  etc.
transitionsLabels = subset(df,Replicate.Name == df$Replicate.Name[1], select=c(8,7,1)) %>% apply(.,1,paste0,collapse=":") %>%  as.vector

# create feature data  "peptide"          "proteinName"      "Precursor.Mz"     "Precursor.Charge" "Product.Mz"       "Product.Charge"   "Fragment.Ion"     "sample"           "isFiltered"  
transFetureData = subset(df,Replicate.Name == df$Replicate.Name[1])
transFetureData =  subset(transFetureData, select=which(!(names(transFetureData) %in%  c("Replicate.Name","Area","Background","Peak.Rank","Retention.Time"))) )
transFetureData$isFiltered = F
names(transFetureData)[names(transFetureData) == "Peptide.Sequence"] = "peptide"
names(transFetureData)[names(transFetureData) == "Protein.Name"] = "proteinName"
row.names(transFetureData) =transitionsLabels 

### create expression matrix ( each column is a run, i.e. dilution replicate)
transExpressionMatrix = data.frame( row.names=transitionsLabels )
for(repName in unique(df$Replicate.Name)){
  #transExpressionMatrix = cbind(transExpressionMatrix,subset(df,Replicate.Name == repName)$Area )
  transExpressionMatrix = cbind(transExpressionMatrix, df[df$Replicate.Name == repName,]$Area %>% as.character %>% as.numeric )
}
names(transExpressionMatrix) = unique(df$Replicate.Name)
transExpressionMatrix = as.matrix(transExpressionMatrix) 

### create expe design
expDesign = data.frame(row.names=colnames(transExpressionMatrix), isControl=colnames(transExpressionMatrix) %>% grepl("HEK",.), condition=  colnames(transExpressionMatrix) %>% substr(.,1,3)   )
# add concentration per sample to pheno data
expDesign$concentration = c( rep(10^(1:7),2) %>% sort, rep(0,4))
transEset = createExpressionDataset(expressionMatrix=transExpressionMatrix,expDesign=expDesign,featureAnnotations=transFetureData)

### create peptide  eset
peptideEset = rollUp(transEset, featureDataColumnName = "peptide")


# > pData(peptideEset)
# isControl condition concentration
# 002_DC081116             FALSE       002         1e+01
# 002_DC081116_Rep         FALSE       002         1e+01
# 005_DC081116             FALSE       005         1e+02
# 005_DC081116_Rep         FALSE       005         1e+02
# 008_DC081116_C16         FALSE       008         1e+03
# 008_DC081116_Rep_C16     FALSE       008         1e+03
# 011_DC081116_C16         FALSE       011         1e+04
# 011_DC081116_Rep_C16     FALSE       011         1e+04
# 014_DC081116_C16         FALSE       014         1e+05
# 014_DC081116_Rep_C16     FALSE       014         1e+05
# 017_DC081116_C16         FALSE       017         1e+06
# 017_DC081116_Rep_C16     FALSE       017         1e+06
# 020_DC081116_C16         FALSE       020         1e+07
# 020_DC081116_Rep_C16     FALSE       020         1e+07
# HEK05_DC081116            TRUE       HEK         0e+00
# HEK06_DC081116            TRUE       HEK         0e+00
# HEK07_DC081116            TRUE       HEK         0e+00
# HEK08_DC081116            TRUE       HEK         0e+00

### INIT END

testGetLOD <- function(){
	
	cat(" --- testGetLOD --- \n")
  df = data.frame(concentration= pData(peptideEset)$concentration,intensity=exprs(peptideEset)[2,])
  method= "blank"
  stopifnot(round(getLOD(df,method="blank")) == 1928)
  stopifnot(round(getLOD(df,method="low")) == 4012)
	cat(" --- testGetLOD: PASS ALL TEST  --- \n")
	
}

testGetDotProduc = function(){
  
  cat(" --- testGetDotProduc --- \n")

  A = 1:3 
  B = 2*A 
  stopifnot(dotProduct(A,B,norm=T) == 1)
  
  cat(" --- testGetDotProduc: PASS ALL TEST  --- \n")
  
}

testGetAllDotProduc = function(){
  
  cat(" --- testGetAllDotProduc --- \n")
  stopifnot(which.max(getAllDotProduct(transEset)[1,]) == 14)
  cat(" --- testGetAllDotProduc: PASS ALL TEST  --- \n")
  
}


### RUN TESTS

testGetLOD()
testGetDotProduc()
testGetAllDotProduc()








