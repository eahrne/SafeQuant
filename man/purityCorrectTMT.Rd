% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TMT.R
\name{purityCorrectTMT}
\alias{purityCorrectTMT}
\title{Correct channel intensities based on Reporter ion Isotopic Distributions}
\usage{
purityCorrectTMT(tmtData, impurityMatrix = impurityMatrix)
}
\arguments{
\item{tmtData}{data.frame containing tmt channel intensities}

\item{impurityMatrix}{correction matrix}
}
\value{
data.frame of corrected tmt intensities
}
\description{
Correct channel intensities based on Reporter ion Isotopic Distributions
}
\details{
Same method as MSnbase, and described in Breitwieser et al. 2012 (Book Chapter)
}
\note{
No note
}
\examples{
print("No examples")
}
\references{
NA
}
