--- Introduction

SafeQuant allows for statistical validation of differentially expressed (DE) proteins/peptides across two or more sample conditions. The analysis requires LC-MS runs from biological triplicate samples, per condition, aligned using Progenesis LC-MS (NonLinear Dynamics).
For each peptide/protein (identified using a database search tool such as Mascot, Sequest or X!Tandem) an expression ratio is calculated based on the median feature intensity of the biological triplicate samples. Each ratio is subjected to statistical validation accounting for the variation of peptide intensities per condition, as well as the total number of detected features (adjustment for multiple testing). In addition to a simple spreadsheet export the SafeQuant analysis output provide numerous graphics to visualize global differences between sample conditions.

- Getting started
-Download the latest version of SafeQuant.
-Install R, 'The R Project for Statistical Computing'
-Install required R libraries

limma
affy
gplots
optparse

--- Using the tool
SafeQuant is executed from the command line and consists of two separate modules designed to post-process (1) peptide level Progenesis exports (peptideAnalysis.R) and (2) protein level Progenesis exports. Both tools have very similar functionality, where the the peptide level analysis include some additional features with respect to the protein level analysis.
Below follows a brief guide on how to perform a peptide level analysis of differentially expressed proteins. A couple of Test datasets are available in the ./SafeQuant/PeptidesSQAnalysis/tesData/ directory.

- 1) Complete an Progenesis LC-MS analysis workflow. Generate spreadsheet (.csv) of Progenesis results. File -> Export Feature Data .. -> (keep default check-box selection (Imoprtant!)
(When performing a protein level analysis, File -> Export Protein Measurements .. (Default Selection!) )

- 2) Open shell and navigate to the SafeQuant installation directory.

erikahrne$ cd /Users/user/mytools/SafeQuant/

- 3) Execute analysis with default options

erikahrne$ Rscript ./PeptidesSQAnalysis/peptideAnalysis.R -i PeptidesSQAnalysis/testData/peptides1.csv -o PeptidesSQAnalysis/out/

-i path to input file -o path to output directory -l output file label

Note that Windows users may have to provide the full path to Rscript
e.g.
c:\Program Files\R\R-2.14.0\bin\Rscript ./PeptidesSQAnalysis/peptideAnalysis.R -i PeptidesSQAnalysis/testData/peptides1.csv -o PeptidesSQAnalysis/out/

- 4) Display options

erikahrne$ Rscript ./PeptidesSQAnalysis/peptideAnalysis.R -h

- 5) Further explanation of some of the user options:

Data filtering The user has several options for data filtering:
I) Peptide level FDR cutoff. Provided that the Progenesis .csv export file includes decoy entries, spurious peptide identification can be filtered out based on a peptide level FDR cutoff specified by the user (-f option default 0.01)
II) Protein Accession filter The user can choose to focus the DE analysis on peptides belonging to a certain group of Proteins or carrying a defined Post-Translational Modification (e.g. Phosphorylation) (-a and Ðt options)
III) Precursor Mass Error Filter (peptide analysis specific option) The user can choose to exclude all peptide identifications with a precursor mass error exceeding a set threshold (-m option, default 10 ppm)
Data normalization After the initial data filtering steps and replacement of missing intensity values the data is normalized i.e. the total intensity sum of all peptide features is scaled to be equal for all samples.
Missing peptide intensity values can be replaced by a user-defined intensity (-n option default 100).
Analysis of Differential Expression Feature intensities are summed per unique peptide sequence. Modified and non-modified variants of the same peptide are treated as separate peptide entries.
Peptide expression ratios are calculated as the intensity ratios of triplicate medians. Each expression ratio is associated with a *q-value (p-value adjusted for multiple testing)
*a modified t-statistic is calculated using the empirical Bayes method, Smyth (2004), and subsequently adjusted for multiple testing using the Benjamin & Hochberg (1995) method.

- TIP !!!
Processing large peptide level input files can be slow, it may be meaningful to filter out features not annotated with a peptides sequence, prior to SafeQuant analysis. The perl script filterLargeProgenesisPeptideFile.pl ( available in the ./SafeQuant/script/ directory) can be used for this purpose.
erikahrne$ perl filterLargeProgenesisPeptideFile.pl largeFile.csv

---  Output files

-  Graphics
The .pdf file contains various graphics visualizing global differences between sample conditions.

-  Spread sheet (.csv)
The exported .csv file contains all data underlying the plots listed in the .pdf file e.g. median ratios and q-values per peptide/protein and condition comparison. Data columns are tab-delimitated and can be explored using e.g. Microsoft Excel or Open Office Calc.

-  .RData file
Import data structures created by SafeQuant into R.

R> load("C:/Users/data/SafeQuant_Analysis.RData")

Contains the following objects:

PROTEINS
GROUPEDRAWINTDATA
GROUPEDSPECCOUNTDATA
GROUPEDNORMINTDATA
MEDIANNORMINTDATA
RATIOSMEDIANNORMINTPERCOND
QVALUESPERCOND
CVSPERCOND
CONTROLCONDITION
EXPDESIGN

--- Publications

* Large-scale quantitative assessment of different in-solution protein digestion protocols reveals superior cleavage efficiency of tandem Lys-C/trypsin proteolysis over trypsin digestion
Timo Glatter, Christina Ludwig, Erik Ahrne, Ruedi Aebersold, Albert J.R. Heck, and Alexander Schmidt
Journal of Proteome Research Just Accepted Manuscript
