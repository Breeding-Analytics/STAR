\name{ClusterKmeans}
\alias{ClusterKmeans}
\alias{ClusterKmeans.default}
\title{
K-means Cluster Analysis
}
\description{
This function performs the K-means Cluster Analysis for numeric variables only.
}
\usage{
ClusterKmeans(data, var, idVar = NULL, clusterMem = TRUE, descriptiveStatRaw = TRUE, corMatx = TRUE, 
              descriptiveStat = TRUE, kgraph = TRUE, clusterNum = 2, saveMem = TRUE, outputPath = NULL)
}
\arguments{
  	\item{data}{a data frame}
	\item{var}{a non-empty character vector of the variables to be used in the analysis}
	\item{idVar}{a grouping variable}
	\item{clusterMem}{a logical value indicating whether the cluster membership will be displayed.}
  	\item{descriptiveStatRaw}{a}
  	\item{corMatx}{a}
	\item{descriptiveStat}{a logical value indicating whether the ddescriptive statistics will be displayed.}
	\item{kgraph}{a logical value indicating whether the graph will be displayed.}
	\item{clusterNum}{a numeric value which specifies the desired number of groups or clusters. The default value is 2.}
	\item{saveMem}{an option whether to save cluster membership to a file or not.}
  	\item{outputPath}{a path where the output file (graph, text file, csv) will be saved.}
}
\author{
Christoffer Edd Relente
}
\keyword{ multivariate }
