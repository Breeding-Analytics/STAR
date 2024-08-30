\name{MergeData}
\alias{MergeData}
\alias{MergeData.default}
\title{Merge Dataset}
\description{Merge two data frame}
\usage{
MergeData(MasterData, TransactionData, byMaster, byTransact, MasterVarInclude = names(MasterData), TransactVarInclude = names(TransactionData), allMaster = FALSE, allTransact = FALSE)
}
\arguments{
  \item{MasterData}{data frame}
  \item{TransactionData}{data frame}
  \item{byMaster}{character string or vector}
  \item{byTransact}{character string or vector}
  \item{MasterVarInclude}{character string or vector}
  \item{TransactVarInclude}{character string or vector}
  \item{allMaster}{logical}
  \item{allTransact}{logical}
}
\author{AAGulles}
\keyword{dataset}
