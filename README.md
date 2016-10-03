The SafeQuant Package includes methods for analysis of quantitative (LFQ,TMT,HRM) Proteomics data.

### Installation

#### 1) Install Dependencies

A) Install CRAN library dependencies (open R)

	R> install.packages(c("seqinr","gplots","corrplot","optparse","data.table","epiR"))

B) Install BioConductor library dependencies (open R)

	R> source("http://bioconductor.org/biocLite.R")
	R> biocLite(c("limma","affy"))

#### 2) Install SafeQuant from sources

**Option 1, install "master branch" using "devtools"**

Make sure you have a working development environment.

**Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

**Mac**: Install Xcode from the Mac App Store.

**Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

    R> install.packages("devtools")
    R> library(devtools)
    R> install_github("eahrne/SafeQuant")
    
**Option 2, install latest [CRAN](https://CRAN.R-project.org/package=SafeQuant) version**

	R> install.packages("SafeQuant")

#### 3) To run safeQuant.R (Post-process Progenesis LFQ datasets or Scaffold TMT datasets)

A) locate file safeQuant.R (C:\Users\ahrnee-adm\Downloads\SafeQuant\exec\safeQuant.R ) 
This is the SafeQuant main script. Copy it to an appropriate directory, e.g. c:\Program Files\SafeQuant\
	
B) open terminal
To display help options

	> Rscript "c:\Program Files\SafeQuant\safeQuant.R" -h
To run (with minimal arguments)

	> Rscript "c:\Program Files\SafeQuant\safeQuant.R" -i "c:\Program Files\SafeQuant\testData\peptide_measurement.csv" -o "c:\Program Files\SafeQuant\out"

#### Tips

I) If using Progenesis QI we advice running SafeQuant on "Peptide Measurement" Exports. 
- File -> Export Peptide Measurements.  This option is available once you have reached the "Resolve Conflicts" Step in Progenesis QI
- When choosing properties to be included in the exported file check the "Grouped accessions (for this sequence)" check box.

II) When working with Progenesis "Feature Exports" it is advisable to discard all features (rows) not annotated with a peptide, to speed up SafeQuant analysis.
This can be done using the "filterLargeProgenesisPeptideFile.pl" perl script. (C:\Users\ahrnee-adm\Downloads\SafeQuant\exec\filterLargeProgenesisPeptideFile.pl) 

A) install perl (or activePerl for windows http://www.activestate.com/activeperl)
	
B) open terminal

	> perl "C:\Program Files\SafeQuant\filterLargeProgenesisPeptideFile.pl" "C:\Program Files\SafeQuant\testData\features.csv"
This will create a new versions of the feature file called with the extension "_FILTERED" features.csv -> features_FILTERED.csv

#### Basic functionality of the safeQuant.R script

1. **Data Normalization**
	* LFQ
		* Global data normalization by equalizing the total MS1 peak areas  across all LC/MS runs.
	* Isobaric Labeling experiments (TMT or iTRAQ)
		* Global data normalization by equalizing the total reporter ion intensities across all reporter ion channels.
2. **Ratio Calculation**
	* LFQ
		* Summation of MS1 peak areas per peptide/protein and LC-MS/MS run, followed by calculation of peptide/protein abundance ratios. 
	* Isobaric Labeling experiments (TMT or iTRAQ)
		* Summation of reporter ion intensities per peptide/protein and LC-MS/MS run, followed by calculation of peptide/protein abundance ratios. 
3. **Statistical testing for differential abundances**
	* The summarized peptide/protein expression values are used for statistical testing of between condition differentially abundant peptides/proteins. Here, empirical Bayes moderated t-tests is applied, as implemented in the R/Bioconductor limma package (Smyth, 2004). The resulting per protein and condition comparison p-values are subsequently adjusted for multiple testing using the Benjamini-Hochberg method.

Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol, 3 SP -, Article3. http://www.ncbi.nlm.nih.gov/pubmed/16646809

#### Use Case Manual

https://raw.githubusercontent.com/eahrne/SafeQuant/master/inst/manuals/SafeQuant_UseCases.txt

#### .tsv export help

https://github.com/eahrne/SafeQuant/blob/master/inst/manuals/tsv_spreadsheet_help.pdf

#### Package Documentation

https://github.com/eahrne/SafeQuant/blob/master/inst/manuals/SafeQuant-man.pdf

#### Publications

* Ahrne, E. et al. Evaluation and Improvement of Quantification Accuracy in Isobaric Mass Tag-Based Protein Quantification Experiments. J Proteome Res 15, 2537–2547 (2016). https://www.ncbi.nlm.nih.gov/pubmed/27345528 
* Ahrne, E., Molzahn, L., Glatter, T., & Schmidt, A. (2013). Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Proteomics. Journal of Proteome Research Just Accepted Manuscript https://www.ncbi.nlm.nih.gov/pubmed/23794183
* Glatter, T., Ludwig, C., Ahrne, E., Aebersold, R., Heck, A. J. R., & Schmidt, A. (2012). Large-scale quantitative assessment of different in-solution protein digestion protocols reveals superior cleavage efficiency of tandem Lys-C/trypsin proteolysis over trypsin digestion. https://www.ncbi.nlm.nih.gov/pubmed/23017020

**[questions?](mailto:erik.ahrne@unibas.ch)**