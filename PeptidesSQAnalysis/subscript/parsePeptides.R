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
source(paste(BASEDIR,"PeptidesSQAnalysis/functions/CsvParserPeptides.R", sep=""))

#params 
INTSTARTCOLUMN <- 12

#get experimental design
EXPDESIGNOBJ <- parseExpDesign(userOptions$progenesisFilePath, verbose=verbose)
EXPDESIGN <- EXPDESIGNOBJ$expDesign 
NBSAMPLES <- EXPDESIGNOBJ$nbSamples
ISENOUGHREPLICATES <- EXPDESIGNOBJ$isEnoughReplicates
NBCONDITIONS <- length(names(EXPDESIGN))


#parse progenesis file all data 
PROGOUTDATA <- getSummedUpPeptideProgOutData(parseProgOutData(userOptions$progenesisFilePath,verbose=userOptions$verbose)
		,intStartColumn =INTSTARTCOLUMN
		,nbSamples=NBSAMPLES
		, proteaseTarget=userOptions$proteaseTarget
		,verbose=userOptions$verbose
)

#apply filters
FILTEROBJ <- filterProgOutDataPeptides(PROGOUTDATA
		, selectedProteinName=userOptions$selectedProteinName
    	, selectedModifName=userOptions$selectedModifName
		, fdrCutoff=userOptions$fdrCutoff
		, precursorMassFilter = userOptions$precursorMassFilter
		, verbose=userOptions$verbose)
PROGOUTDATA <- FILTEROBJ$progOutData

### test user specified control condition
if(userOptions$selectedControlCond > NBCONDITIONS){
	print(paste("ERROR: Invalid control condition specified ->",userOptions$selectedControlCond  ))
	quit(status=-1)
}

# create data structures UNGROUPEDRAWINTDATA, GROUPEDRAWINTDATA, 

UNGROUPEDRAWINTDATA <- getUngroupedData(PROGOUTDATA,columnStart=9, nbSamples=NBSAMPLES)
#GROUPEDRAWINTDATA <- getCondGroupedSampleData(UNGROUPEDRAWINTDATA,expDesign=EXPDESIGN)

UNGROUPEDSPECCOUNTDATA <- getUngroupedData(PROGOUTDATA,columnStart=9+NBSAMPLES, nbSamples=NBSAMPLES)

###################################### PARSE CSV FILE END
