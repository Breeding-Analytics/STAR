\name{CompareMeans}
\alias{CompareMeans}
\alias{CompareMeans.default}
\title{t-Test}
\description{descriptionTemp}
\usage{
CompareMeans(data, varX, varY = NULL, grp = NULL, testVal = 0, procedure = c("one", "paired", "independent"), alternative = c("two.sided", "greater", "less"), statistics = FALSE, CI = FALSE, confLevel = 0.95, normality = NULL, alpha = 0.05)
}
\arguments{
  \item{data}{name of data frame}
  \item{varX}{numeric variable}
  \item{varY}{numeric variable}
  \item{grp}{grouping variable}
  \item{testVal}{hypothesized value}
  \item{procedure}{statistical procedure}
  \item{alternative}{alternative hypothesis}
  \item{statistics}{logical}
  \item{CI}{logical}
  \item{confLevel}{level of confidence}
  \item{normality}{NULL or normality test procedure}
  \item{alpha}{level of significance}
}
\author{AAGulles}
\keyword{ means }
