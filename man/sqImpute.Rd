% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExpressionAnalysis.R
\name{sqImpute}
\alias{sqImpute}
\title{Impute missing values}
\usage{
sqImpute(eset, method = "gmin", rowmax = 0.3)
}
\arguments{
\item{eset}{ExpressionSet}

\item{method}{c("knn","ppca","gMin",lMin)}

\item{rowmax}{The maximum percent missing data allowed in any row to apply ppca and knn (if more missing values impute gmin). default 0.3}
}
\value{
ExpressionSet
}
\description{
Impute missing values
}
\details{
\itemize{
\item{gMin: half global minimum (0.1 percentile)}
\item{lMin: half local minimum}
\item{gMean: half global mean}
\item{lMean: half local mean}
\item{knn: Nearest neighbour averaging, as implemented in the impute::impute.knn function}
\item{ppca: An iterative method using a probabilistic model to handle missing values, as implemented in the pcaMethods::pca function.}
}
}
\note{
No note
}
\examples{
print("No examples")
}
\references{
Accounting for the Multiple Natures of Missing Values in Label-Free Quantitative Proteomics Data Sets to Compare Imputation Strategies, Lazar et al (2016), \url{http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00981}
}
\seealso{
No note
}
