#!/usr/bin/Rscript


# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### GET BASEDIR
initial.options <- commandArgs(trailingOnly = FALSE)
BASEDIR <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
BASEDIR <- gsub("ProteinsSQAnalysis","",BASEDIR,)
BASEDIR <- gsub("^\\.$","../",BASEDIR,perl=T)
remove(initial.options)
### GET BASEDIR END

ISPEPTIDEANALYSIS <- F

source(paste(BASEDIR,"scripts/runAnalysis.R", sep=""))