\name{PCA}
\alias{PCA}
\alias{PCA.default}
\title{Principal Component Analysis}
\description{
Performs Principal Component Analysis for numeric variables only.
}
\usage{
PCA(data, var, idVar = NULL, descriptiveStat = TRUE, corMatx = TRUE, covMatx = TRUE, matx = c("corr","cov"), 
    transform = c("zerocenter", "unitvar", "none"), saveScore = TRUE, outputPath = NULL, scatterMatx = TRUE, 
    scree = TRUE, biplot = TRUE, pcaplot = TRUE, useIdVar = FALSE, pChars, pSizes, pCol, showLeg = FALSE, legTitle, 
    legPos = "bottomright", legNcol = 1, axesNum = 3)
}
\arguments{
  	\item{data}{a data frame}
  	\item{var}{a non-empty character vector of the variables to be used in the analysis}
	\item{idVar}{a grouping variable}
	\item{descriptiveStat}{logical}
  	\item{corMatx}{logical}
  	\item{covMatx}{logical}
  	\item{matx}{a character string which specifies whether the matrix to be used is correlation or covariance.}
  	\item{transform}{a character string which specifies whether the data transformation to be used is zero- centered, unit variance or no transformation.}
	\item{saveScore}{an option whether to save pc scores to a file or not.}
  	\item{outputPath}{a path where the output file (graph, text file, csv) will be saved.}
	\item{scatterMatx}{g}
	\item{scree}{g}
	\item{biplot}{g}
	\item{pcaplot}{g}
	\item{useIdVar}{a}
	\item{pChars}{g}
	\item{pSizes}{g}
	\item{pCol}{g}
	\item{showLeg}{g}
	\item{legTitle}{g}
	\item{legPos}{g}
	\item{legNcol}{g}
	\item{axesNum}{g}
}
\details{
	The function performs the analysis using prcomp from package stats. 
	If transform is normal, the data are standardized while if transform is center, 
	the data are zero centered only before the analysis takes place.
}
\author{
Christoffer Edd Relente
}
\keyword{ multivariate }