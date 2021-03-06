% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExpressionAnalysis.R
\name{rollUpDT}
\alias{rollUpDT}
\title{Roll up feature intensites per unique colum combination}
\usage{
rollUpDT(eset, method = "sum", featureDataColumnName = c("proteinName"))
}
\arguments{
\item{eset}{ExpressionSet}

\item{method}{"sum", "mean" or "top3"}

\item{featureDataColumnName}{vector of column names e.g. peptide or proteinName}
}
\value{
ExpressionSet object
}
\description{
Roll up feature intensites per unique colum combination
}
\details{
featureDataColumnName = c("peptide","charge","ptm"), method= c("sum"), sums up intensities per peptie modification charge state
}
\note{
No note
}
\examples{
print("No examples")
}
\references{
No references
}
