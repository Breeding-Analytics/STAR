\name{ClusterDivisive}
\alias{ClusterDivisive}
\alias{ClusterDivisive.default}
\title{
Divisive Hierarchical Cluster Analysis
}
\description{
This function performs the Divisive Hierarchical Cluster Analysis for Numeric variables only.
}
\usage{
ClusterDivisive(data, var = NULL, idVar =NULL, sbinVar = NULL, abinVar = NULL, ofactorVar = NULL, factorVar = NULL, 
                stand = TRUE, distance = c("Euclidean", "Maximum", "Manhattan", "Minkowski", "Canberra", "Binary", 
		"Simple Matching","Sokal & Sneath", "Hamann coefficient", "Jaccard", "Dice", "Gower"), distMatrix = FALSE,
		copDistance = FALSE, clusterMem = TRUE, descriptiveStatRaw = TRUE, corMatx = TRUE, scatterMatx = TRUE, 
		descriptiveStat = TRUE, dendrogram = TRUE, clusterBox = TRUE, clusterNum = 2, saveMem = TRUE, 
                outputPath = NULL)
}
\arguments{
  	\item{data}{a data frame}
  	\item{var}{a non-empty character vector of the variables to be used in the analysis}
	\item{idVar}{a grouping variable}
	\item{sbinVar}{a}
	\item{abinVar}{a}
	\item{ofactorVar}{a}
	\item{factorVar}{a}
	\item{stand}{a logical value indicating whether the variables are to be standardized or not. 
	The default value is TRUE.}
	\item{distance}{a character string which specifies the distance measure to be used.}
	\item{distMatrix}{a logical value indicating whether the distance matrix will be displayed. The default value is FALSE.}
	\item{copDistance}{a}
	\item{clusterMem}{a logical value indicating whether the cluster membership will be displayed. The default value is TRUE.}
	\item{descriptiveStatRaw}{a}
	\item{corMatx}{a}
	\item{scatterMatx}{a}
	\item{descriptiveStat}{a logical value indicating whether the descriptive statistics will be displayed. The default value is TRUE.}
	\item{dendrogram}{a logical value indicating whether the dendrogram will be displayed. The default value is TRUE.}
	\item{clusterBox}{a logical value indicating whether the dendrogram will be highlighted. The default value is TRUE.}
  	\item{clusterNum}{a numeric value which specifies the desired number of groups or clusters. The default value is 2.}
	\item{saveMem}{an option whether to save cluster membership to a file or not.}
  	\item{outputPath}{a path where the output file (graph, text file, csv) will be saved.}
}
\author{Christoffer Edd Relente}
\keyword{ multivariate }

