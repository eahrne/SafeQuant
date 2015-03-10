--- Introduction

The SafeQuant Package includes methods for analysis of Proteomics "Quant data".
More documentation to come.

1) Install Dependencies

	A) Install CRAN library dependencies (open R)
	R> install.packages(c("seqinr","gplots","corrplot","optparse","data.table"))

	B) Install BioConductor library dependencies (open R)
	R> source("http://bioconductor.org/biocLite.R")
	R> biocLite(c("limma","affy"))

2) Install SafeQuant from sources

	A) Download Zip from github
	https://github.com/eahrne/SafeQuant

	B) Unzip SafeQuant-v.X.XX-pack.zip
		
	C) Install SafeQuant (open R)
	Assuming that the SafeQuant source directory is at "C:\\Users\\ahrnee-adm\\Downloads\\"
	The SafeQuant source directory is the directory containing the NAMESPACE and DESRIPTION files etc.
	R> install.packages("C:\\Users\\ahrnee-adm\\Downloads\\SafeQuant-v.X.XX-pack",type="source", repos=NULL)

3) To run safeQuant.R (Post-process Progenesis LFQ datasets or Scaffold TMT datasets)

	A) locate file safeQuant.R (C:\Users\ahrnee-adm\Downloads\SafeQuant-v.2.01-pack\exec\safeQuant.R ) 
	This is the SafeQuant main script. Copy it to an appropriate directory, e.g. c:\Program Files\SafeQuant\
	
	B) open terminal
 	# To display help options
	> Rscript "c:\Program Files\SafeQuant\safeQuant.R" -h
 	# To run (with minimal arguments)
	>Rscript "c:\Program Files\SafeQuant\safeQuant.R" -i "c:\Program Files\SafeQuant\testData\proteins.csv" -o "c:\Program Files\SafeQuant\out"

Tips

I) When working with Progenesis "Feature Exports" it is advisable to discard all features (rows) not annotated with a peptide, to speed up SafeQuant analysis.
This can be done using the "filterLargeProgenesisPeptideFile.pl" perl script. (C:\Users\ahrnee-adm\Downloads\SafeQuant-v.2.01-pack\exec\filterLargeProgenesisPeptideFile.pl) 

	A) install perl (or activePerl for windows http://www.activestate.com/activeperl)
	
	B) open terminal
	> perl "C:\Program Files\SafeQuant\filterLargeProgenesisPeptideFile.pl" "C:\Program Files\SafeQuant\testData\features.csv"
	# this will create a new veriosn of the feature file called with the extension "_FILTERED" features.csv -> features_FILTERED.csv

II) If using Progenesis QI we advice running SafeQuant on "Peptide Measurement" Exports. 
    - File -> Export Peptide Measurements.  This option is available once the you reached the "Resolve Conflicts" Step in Progenesis QI
    - When choosing properties to be included in the exported file check the "All accessions (for this sequence)" check box.

--- Publications

* Large-scale quantitative assessment of different in-solution protein digestion protocols reveals superior cleavage efficiency of tandem Lys-C/trypsin proteolysis over trypsin digestion
Timo Glatter, Christina Ludwig, Erik Ahrne, Ruedi Aebersold, Albert J.R. Heck, and Alexander Schmidt
Journal of Proteome Research Just Accepted Manuscript
