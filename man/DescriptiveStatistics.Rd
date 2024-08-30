\name{DescriptiveStatistics}
\alias{DescriptiveStatistics}
\alias{DescriptiveStatistics.default}
\title{Descriptive Statistics}
\description{
  Displays univariate summary statistics for several variables per levels of a 
  grouping variable in a single table.
}
\usage{
  DescriptiveStatistics(data, var, grp = NULL, statistics = c("nnmiss", "mean", "sd"))
}
\arguments{
  \item{data}{a non-null character string which indicates the name of the data frame 
              containing the variables specified in \code{var} and \code{grp} }
  \item{var}{a non-empty character vector which indicates the name of the variable(s) 
             where summary statistics will be computed}
  \item{grp}{\code{NULL} (default) or a character vector which indicates the name of the
             variable(s) whose value define the subgroup for generating the descriptive statistic}
  \item{statistics}{a non-empty character vector which indicates the statistics to be computed. 
                    Must be at least one of the following: \dQuote{n}, \dQuote{nnmiss}, \dQuote{nmiss},
		    \dQuote{sum}, \dQuote{css}, \dQuote{ucss}, \dQuote{se.skew}, \dQuote{se.kurtosis},
		    \dQuote{range}, \dQuote{iqr}, \dQuote{var}, \dQuote{sd}, \dQuote{se.mean},
		    \dQuote{cv}, \dQuote{mean}, \dQuote{median}, \dQuote{mode}, \dQuote{min},
		    \dQuote{max}, \dQuote{q1}, \dQuote{q3}, \dQuote{skew} and/or \dQuote{kurtosis}.}
}
\author{Alaine A. Gulles}
\keyword{statistics}

