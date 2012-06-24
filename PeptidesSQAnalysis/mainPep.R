#!/usr/bin/Rscript


# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### GET BASEDIR
initial.options <- commandArgs(trailingOnly = FALSE)
BASEDIR <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
BASEDIR <- gsub("PeptidesSQAnalysis","",BASEDIR,)
BASEDIR <- gsub("^\\.$","../",BASEDIR,perl=T)
remove(initial.options)
### GET BASEDIR END

ISPEPTIDEANALYSIS <- T

source(paste(BASEDIR,"scripts/runAnalysis.R", sep=""))