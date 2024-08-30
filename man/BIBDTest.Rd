\name{BIBDTest}
\alias{BIBDTest}
\alias{BIBDTest.default}
\title{Balanced Incomplete Block Design}
\description{
 Perform analysis of variance for balanced incomplete block design where 
 the treatments are treated as fixed factor.
}
\usage{
BIBDTest(data, respvar, trmt, block, set = NULL, method = NULL, descriptive = FALSE, normality = FALSE, homogeneity = FALSE, alpha = 0.05, outputPath = NULL)
}
\arguments{
  \item{data}{a data frame}
  \item{respvar}{a character string denoting the variable name}
  \item{trmt}{a character string denoting the treatment}
  \item{block}{a character string denoting the block}
  \item{set}{a character string denoting the by variable}
  \item{method}{a character string}
  \item{descriptive}{logical; denoting whether the descriptive statistics will be displayed}
  \item{normality}{logical; denoting whether the descriptive statistics will be displayed}
  \item{homogeneity}{logical; denoting whether the descriptive statistics will be displayed}
  \item{alpha}{a numeric value denoting the level of significance}
  \item{outputPath}{path where graphical output will be save}
}
\author{
 AAGulles
}
\keyword{analysis}
