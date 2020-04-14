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

proteinSeq1 <- "MSAGSSCSQTPSRAIPTRRVALGDGVQLPPGDYSTTPGGTLFSTTPGGTRIIYDRKFLMECRNSPVAKTPPKDLPAIPGVTSPTSDEPPMQASQSQLPSSPEDKRAGGEESQFEMDI"
proteinSeq2 <- "MVKKSRRRGAAQWAAVRAQAGLTATDENEDDLGLPPSPGDSSYYQDQVDEFHEARSRAVLAKGWNEVESGEEDGDEEEE"
proteinSeq3 <- "MSERMRPVVVDLPTSASSSMKVNG"
proteinSeq4 <- "RKR"

### read protein db
proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

# read phospho motif example file
phosphoMotifs = read.csv(phosphoMotifFile, sep="\t")
### INIT END

### TEST FUNCTIONS

testIsCon <- function(){
	cat("--- testIsCon: --- \n")
	stopifnot(sum(isCon(fData(eset)$proteinName)) == 0)
	cat("--- testIsCon: PASS ALL TEST --- \n")
}

testIsDecoy <- function(){
	cat("--- testIsDecoy: --- \n")
	stopifnot(sum(isDecoy(fData(eset)$proteinName)) == 198)
	cat("--- testIsDecoy: PASS ALL TEST --- \n")
}

testGetIdQvals <- function(){
	cat("--- testGetIdQvals: --- \n")
	stopifnot(max(getIdLevelQvals(fData(eset)$idScore, isDecoy(fData(eset)$proteinName))) < 1)
	cat("--- testGetIdQvals: PASS ALL TEST --- \n")
}


testAddIdQvalues <- function(){

	cat("--- testAddIdQvalues: --- \n")
	stopifnot(all.equal(round(max(fData(addIdQvalues(eset))$idQValue ),2),0.28))
	stopifnot(all.equal(round(max(fData(addIdQvalues(rollUp(eset, featureDataColumnName = "proteinName"  )))$idQValue),2),0.29))
	cat("--- testAddIdQvalues: PASS ALL TEST --- \n")

	#progenesisPeptideCsvFile3 <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/PeptidesSQAnalysis/peptides5.csv"
	#d <- parseProgenesisPeptideCsv(file=progenesisPeptideCsvFile3,expDesign=getExpDesignProgenesisCsv(progenesisPeptideCsvFile3))
	#e <- rollUp(d, method="sum", isProgressBar=T, featureDataColumnName= c("proteinName"))
	#x <- addIdQvalues(e)
	#
	#fData(x)$idQValue
	#isPtm <- nchar(as.character(fData(x)$ptm)) > 0
	#plot(sort(fData(x)$idScore[isPtm]),sort(fData(x)$idQValue[isPtm],decreasing=T), type="l" )
	#lines(sort(fData(x)$idScore[!isPtm]),sort(fData(x)$idQValue[!isPtm],decreasing=T), col=2)
	#
	#
	#d <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisProteinCsvFile1))
	#x <- addIdQvalues(d)
	#fData(x)$idQValue
	#plot(sort(fData(x)$idScore),sort(fData(x)$idQValue,decreasing=T), type="l" )
	#
	#hist(fData(x)$idQValue)


}

testGetScoreCutOff <- function(){
	cat("--- testGetScoreCutOff: --- \n")
	stopifnot(all.equal(round(getScoreCutOff(fData(eset)$idScore, isDecoy(fData(eset)$proteinName)),2), 12.07))
	cat("--- testGetScoreCutOff: PASS ALL TEST --- \n")
}

testGetModifProteinCoordinates <- function(){

	ptm <- "[N-term] Acetyl (Protein N-term)|[6] Oxidation (M)"
	peptide <- "SSDAEMAVFGEAAPYLR"
	protein <-  "MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK"

	#debug(getModifProteinCoordinates)

	cat("--- testGetModifProteinCoordinates: --- \n")
	mc <- getModifProteinCoordinates(ptm,peptide,protein)
	stopifnot(9 == mc[1])
	stopifnot(14 == mc[2])
	stopifnot(all.equal(c(34,49),getModifProteinCoordinates(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1)))
	stopifnot(37 == getModifProteinCoordinates(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2))

	#### SCAFFOLD PTM FORMAT
	proteinSeq <-  "MMMMMMMMMMETPSPRPPPMRHRSSRSP"
	peptideSeq <- "ETPSPRPPPMR"
	modifAnnot <- "T2 Phospho, S4 Phospho, M10 Oxidation"
	stopifnot(all.equal(c(12,14,20),  getModifProteinCoordinates(modifAnnot,peptideSeq,proteinSeq, format=2)))


	cat("--- testGetModifProteinCoordinates: PASS ALL TEST --- \n")
}

testGetMotifX <- function(){

	cat("--- testGetMotifX: --- \n")

	modifPos1 <- getModifProteinCoordinates(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1)
	modifPos2 <- getModifProteinCoordinates(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2)

	stopifnot("PGDYS*TTPG" == getMotifX(modifPos1,fData(eset)$peptide[1],proteinSeq1,4)[1])
	stopifnot("GLPPS*PGDS" == getMotifX(modifPos2,fData(eset)$peptide[2],proteinSeq2,4))

	#stopifnot("PGDYS*TTPG" == getMotifX(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1,4)[1])
	#stopifnot("GLPPS*PGDS" == getMotifX(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2,4))

	cat("--- testGetMotifX: PASS ALL TEST --- \n")

}

testAddPTMCoord <- function(){

	cat("--- testAddPTMCoord: --- \n")
	eset <- .addPTMCoord(eset,proteinDB,motifLength = 4)
	stopifnot( "PGDYS*TTPG,TPGGT*RIIY"  ==  fData(eset)$motifX[1] )
	cat("--- testAddPTMCoord: PASS ALL TEST --- \n")

}


testSetFilter <- function(){

	### peptide data example

	idQValueThrs <- 0.01
	pMassTol <- c(-1,1)
	ptmRegExp <- "?"
	proteinNameRegExp <- "?"

	filter <- data.frame(

			fData(addIdQvalues(eset))$idQValue >= idQValueThrs # id score
			,(fData(eset)$pMassError < pMassTol[1]) | (fData(eset)$pMassError > pMassTol[2]) # precursor mass tolerance
			,isDecoy(fData(eset)$proteinName)	# decoy
			,isCon(fData(eset)$proteinName)	# contaminants
			,!(regexpr(proteinNameRegExp,fData(eset)$proteinName,ignore.case=T) == 1) # protein ac
			,!(regexpr(ptmRegExp,fData((eset))$ptm,ignore.case=T) == 1 ) # ptm

	)

	cat("--- testSetFilter: --- \n")

	stopifnot(sum(fData(.setFilter(eset,filter=filter))$isFiltered) == 230)
	stopifnot(sum(fData(.setFilter(eset,filter=cbind(filter,rep(T,nrow(eset)))))$isFiltered) == 900)

	cat("--- testSetFilter: PASS ALL TEST --- \n")

}

testGetPeptides <- function(){

	stopifnot(length(getPeptides(proteinSeq3)) == 3)
	stopifnot(paste(getPeptides(proteinSeq3),collapse="") == proteinSeq3)
	stopifnot(length(getPeptides(proteinSeq3,proteaseRegExp=.getProteaseRegExp("lys-c"))) == 2)

	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=0)) == 2)
	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=1)) == 4)
	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=2)) == 5)

	### digest whole db
	if(F){

		digestedDB <- list()

		for(i in 1:length(names(proteinDB))){
		#for(i in 1:100){

			cat(i,"/",length(proteinDB),"\n")
			ac <- names(proteinDB)[i]

			for(peptide in getPeptides(proteinDB[[ac]])){

				digestedDB[[peptide]] <- c(digestedDB[[peptide]],ac)

			}

		}
	}
}

testGetNbDetectablePeptides <- function(){

	cat("--- testGetNbDetectablePeptides: --- \n")
	stopifnot(getNbDetectablePeptides(getPeptides(proteinSeq1)) == 8)
	stopifnot(getNbDetectablePeptides(getPeptides(proteinSeq2),peptideLength=c(-Inf,Inf)) == length(getPeptides(proteinSeq2)))
	cat("--- testGetNbDetectablePeptides: PASS ALL TEST --- \n")

}

testGetNbMisCleavages <- function(){

	cat("--- testGetNbMisCleavages: --- \n")
	peptide <- c("PEPTIDEK","PERTIDEK","PERKTIDEK","PERTIDE","RRPERPKK")
	protease <- "trypsin"
	stopifnot( all.equal(c(0,1,2,1,2),getNbMisCleavages(peptide)))
	cat("--- testGetNbMisCleavages:  PASS ALL TEST --- \n")

}

testGetPeptidePerProtein <- function(){

	cat("--- testGetPeptidePerProtein: --- \n")
	stopifnot(getNbPeptidesPerProtein(eset)["prot52"] == 3)
	cat("--- testGetPeptidePerProtein: PASS ALL TEST --- \n")

}

testSetPeptidePerProtein <- function(){

	cat("--- testSetPeptidePerProtein: --- \n")
	stopifnot("nbPeptides" %in% names(fData(setNbPeptidesPerProtein(eset))))
	cat("--- testSetPeptidePerProtein: PASS ALL TEST --- \n")

}

testGetMeanCenteredRange <- function(){

	cat("--- testGetMeanCenteredRange: --- \n")
	stopifnot( round(mean(getMeanCenteredRange(fData(eset)$pMassError)),3) == round(mean(fData(eset)$pMassError),3) )
	cat("--- testGetMeanCenteredRange: PASS ALL TEST --- \n")

}

testIsStrippedACs <- function(){

	cat("--- testIsStrippedACs: --- \n")
	stopifnot(!isStrippedACs( sample(names(proteinDB),100) ) )
	stopifnot(isStrippedACs(fData(eset)$proteinName))
	cat("--- testIsStrippedACs: PASS ALL TEST --- \n")
}

testStripACs <- function(){

	cat("--- testStripACs: --- \n")
	stopifnot(isStrippedACs(stripACs(sample(names(proteinDB),100))))
	cat("--- testStripACs: PASS ALL TEST --- \n")
}

testGetAAProteinCoordinates <- function(){

  cat("--- testGetAAProteinCoordinates: --- \n")
  peptide <- "SSDAEMAVFGEAAPYLR"
  protein <-  "MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK"
  stopifnot(length(getAAProteinCoordinates("SSDAEMAVFGEAAPYLR","MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK","S")) == 2)
  stopifnot(getAAProteinCoordinates("SSDAEMAVFGEAAPYLR","MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK","Y") == 23)
  cat("--- testGetAAProteinCoordinates: PASS ALL TEST --- \n")

}

testGetMotifFreq = function(){

  # VEVNTNSGEIIHK -> VEVNTNpSGEIIHK
  motifs = gsub("(.{6})([STY].{6})",paste0("\\1","p","\\2"),phosphoMotifs$motif)
  cat("--- testGetMotifFreq: --- \n")
  motifFreq = getMotifFreq( motifs )
  stopifnot(sum(motifFreq$nbMatchesPerMotif == 0) == 170)
  cat("--- testGetMotifFreq: PASS ALL TEST --- \n")

  # rownames(motifFreq) = motifFreq$motif
  #
  # tmp = table(motifFreq$motif)[table(motifFreq$motif) > 1]
  #
  # unique(kinaseMotif[ kinaseMotif$motif %in% names(tmp),])

}


testGetKinaseFreq = function(){

  # VEVNTNSGEIIHK -> VEVNTNpSGEIIHK
  motifs = gsub("(.{6})([STY].{6})",paste0("\\1","p","\\2"),phosphoMotifs$motif)
  cat("--- testGetKinaseFreq: --- \n")
  kinaseStats = getKinaseFreq( motifs[1:10] )
  stopifnot(sum(kinaseStats) == 45)
  cat("--- testGetKinaseFreq: PASS ALL TEST --- \n")

}

testGetKinases = function(){

  # VEVNTNSGEIIHK -> VEVNTNpSGEIIHK
  motifs = gsub("(.{6})([STY].{6})",paste0("\\1","p","\\2"),phosphoMotifs$motif)
  cat("--- testGetKinases: --- \n")
  stopifnot(    nrow(getKinases( motifs[1] )) == 14)
  cat("--- testGetKinases: PASS ALL TEST --- \n")

}

testGetGeneName = function(){
  ### parse gene name
  description1 = "ATP synthase subunit beta OS=Salmonella typhimurium (strain SL1344) GN=atpD"
  description2 = "14-3-3 protein zeta/delta OS=Mus musculus GN=Ywhaz PE=1 SV=1"
  description3 = "14-3-3 protein zeta/delta OS=Mus musculus PE=1 SV=1"
  description4 = ""
  description5 = NA

  cat("--- testGetGeneName: --- \n")

  all = c(description1,description2,description3,description4,description5)
  res = getGeneName(all)
  stopifnot( all(na.omit(res) == c("atpD","Ywhaz")))

  cat("--- testGetGeneName: PASS ALL TEST --- \n")

}


testGetAccessionNumber = function(){

  cat("--- testAccessionNumber: --- \n")

  proteinName = c("sp|A0JLT2|MED19_HUMAN","sp|A0MZ66-1|SHOT1_HUMAN", "myProtein", "P25665","tr|B4DRR0|B4DRR0_CON-HUMAN","tr|A0A022YWF9|E1WHQ6_SALTS","tr|E1WHQ6|E1WHQ6_SALTS,sp|A0A022YWF9|MED19_HUMAN",NA,"cust|Tet1CD|AlainWeber","P0CG47;P0CG48;P62979;P62987")
#
#   stripped = str_replace_all(proteinName, "(^.{2}\\|)|(\\|.*)","")
#   # make sure UniProtKB accession numbers consist of 6 or 10 alphanumerical characters
#   str_extract(stripped, "^[A-Z][0-9][A-Z 0-9]{3}[0-9][A-Z 0-9]{0,4}(\\-){0,1}[0-9]{0,2}$")
#
  res = getAccessionNumber(proteinName)
  stopifnot( all( c("A0JLT2","A0MZ66-1",NA,"P25665","B4DRR0","A0A022YWF9","E1WHQ6",NA,NA,"P0CG47")  ==  res) %>% na.omit )

  cat("--- testAccessionNumber: PASS ALL TEST --- \n")

}

testGetGoTermDF = function(){

	cat("--- testGetGoTermDF: --- \n")

	taxId = 83333 # ecoli
	acs = c("P00350","P00363","P00370")

	goTermDf = getGoTermDF(taxId,acs)

	stopifnot(nrow(goTermDf) == 3)

	cat("--- testGetGoTermDF: PASS ALL TEST --- \n")

}





### TEST FUNCTIONS END

### TESTS

testIsCon()
testIsDecoy()
testGetIdQvals()
testAddIdQvalues()
testGetScoreCutOff()
testGetModifProteinCoordinates()
testGetMotifX()
testAddPTMCoord()
testSetFilter()
testGetPeptides()
testGetNbDetectablePeptides()
testGetNbMisCleavages()
testGetPeptidePerProtein()
testSetPeptidePerProtein()
testGetMeanCenteredRange()
testIsStrippedACs()
testStripACs()
testGetAAProteinCoordinates()
testGetMotifFreq()
testGetKinaseFreq()
testGetKinases()
testGetGeneName()
testGetAccessionNumber()
testGetGoTermDF()
### TESTS END


##### phospho kinase motif analysis

