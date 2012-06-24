# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### CREATE DATA STRUCTURES
#PROGOUTDATA
#EXPDESIGN
#NBSAMPLES
#ISENOUGHREPLICATES
#NBCONDITIONS
#UNGROUPEDRAWINTDATA
#UNGROUPEDSPECCOUNTDATA

######################################  PARSE CSV FILE
source(paste(BASEDIR,"ProteinsSQAnalysis/functions/CsvParserProteins.R", sep=""))

#params 
INTSTARTCOLUMN <- 10

#parse progenesis file all data 
PROGOUTDATA <- parseProgOutData(userOptions$progenesisFilePath,verbose=userOptions$verbose)

#apply filters
FILTEROBJ <- filterProgOutDataProteins(PROGOUTDATA,selectedProteinName=userOptions$selectedProteinName 
		,minNbPeptidesPerProt=userOptions$minNbPeptidesPerProt ,fdrCutoff=userOptions$fdrCutoff, verbose=userOptions$verbose)
PROGOUTDATA <- FILTEROBJ$progOutData

#get experimental design
EXPDESIGNOBJ <- parseExpDesign(userOptions$progenesisFilePath, verbose=verbose)
EXPDESIGN <- EXPDESIGNOBJ$expDesign 
NBSAMPLES <- EXPDESIGNOBJ$nbSamples
ISENOUGHREPLICATES <- EXPDESIGNOBJ$isEnoughReplicates
NBCONDITIONS <- length(names(EXPDESIGN))

### test user specified control condition
if(userOptions$selectedControlCond > NBCONDITIONS){
	print(paste("ERROR: Invalid control condition specified ->",userOptions$selectedControlCond  ))
	quit(status=-1)
}

# create data structures UNGROUPEDRAWINTDATA, GROUPEDRAWINTDATA, 
UNGROUPEDRAWINTDATA <- getUngroupedData(PROGOUTDATA,columnStart=INTSTARTCOLUMN+NBSAMPLES, nbSamples=NBSAMPLES)
GROUPEDRAWINTDATA <- getCondGroupedSampleData(UNGROUPEDRAWINTDATA,expDesign=EXPDESIGN)

UNGROUPEDSPECCOUNTDATA <- getUngroupedData(PROGOUTDATA,columnStart=INTSTARTCOLUMN+2*NBSAMPLES, nbSamples=NBSAMPLES)

###################################### PARSE CSV FILE END