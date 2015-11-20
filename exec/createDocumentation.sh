#!/bin/sh
Rscript /Users/ahrnee-adm/dev/R/workspace/SafeQuant/exec/roxygenize.R
R CMD Rd2pdf ../../SafeQuant --output=/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/manuals/SafeQuant-man.pdf --force --no-preview