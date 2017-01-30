#'   Kinase motifs
#'
#'  Human Protein Reference Database Serine/Threonine motifs
#'  http://www.hprd.org/serine_motifs
#'  The variables are as follows:
#'
#' @format A data frame with 175 rows and 2 variables:
#' \describe{
#'   \item{motif}{kinase motif)}
#'   \item{kinase}{kinase}
#' }
"kinaseMotif"

 kinaseMotif =  read.csv(file="~/dev/R/workspace/CristinaPhospho/data/motifsEtc/HPRD_kinas_motifs.csv")
 # discard double phoshpho e.g. pS*PXX[pS/pT]
 selDoubleS = lapply(kinaseMotif$motif,  function(t){  gregexpr("pS",t)} %>%  unlist %>% length   ) %>% unlist %>% as.vector > 1
 selDoubleT = lapply(kinaseMotif$motif,  function(t){  gregexpr("pT",t)} %>%  unlist %>% length   ) %>% unlist %>% as.vector > 1
 kinaseMotif = kinaseMotif[!(selDoubleS | selDoubleT | grepl("\\*",kinaseMotif$motif) ),]
 kinaseMotif$regExpr = kinaseMotif$motif %>% gsub(("X"),".",.) %>% gsub("(p)([STY])","(p\\2)",.) %>%  gsub("\\)\\/\\(",")|(",.) %>%  gsub("\\[\\(","((",.)  %>%  gsub("\\)\\]","))",.) %>% gsub("\\/","",.)
 kinaseMotif$kinase = gsub("Kinase","kinase",kinaseMotif$kinase)
 kinaseMotif = unique(kinaseMotif)
 save(kinaseMotif, file="~/dev/R/workspace/SafeQuant/data/kinaseMotif.rda", row.names=F)
# grepl("(pS)|(pT)","XXXpSXXX", ignore.case = FALSE)
# grepl("KKKKKK((pS)|(pT)).","KKKKKKpTX")
# grepl("R.R..((pS)|(pT))[FL]","RXRXXpSL")
# grepl("R.R..((pS)|(pT))[FL]X","RXRXXpSL")
# grepl("R.R..((pS)|(pT))[FL]","RXRXXpTLX")
# grepl("...","KK")