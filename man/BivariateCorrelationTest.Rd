\name{BivariateCorrelationTest}
\alias{BivariateCorrelationTest}
\alias{BivariateCorrelationTest.default}
\title{Test of Association between two variables}
\description{
  Perform test of association between two variables using the Pearson's product moment
  correlation coefficient, Kendall's tau and Spearmans's rho.
}
\usage{
  BivariateCorrelationTest(data, var, method = "pearson", alternative = "two.sided", statistics = FALSE)
}
\arguments{
  \item{data}{a non-null character string which indicates the name of the data frame containing
              the variables specified in \code{var}}
  \item{var}{a non-empty character vector which indicates the name of the variable(s) where test of
             association will be performed}
  \item{method}{a character vector indicating which correlation coefficient is to be used for the test
		and must be at least one of the following: "pearson" (default), "kendall" and/or "spearman"}
  \item{alternative}{a character string indicating the alternative hypothesis to be used and must be one
                     of the following: "two.sided" (default), "greater" or "less"}
  \item{statistics}{logical indicating whether summary statistics will be computed for each of the variables
		    in the object \code{var}}
}
\value{
  A list containing the following components:
  \item{procedure }{a character string indicating the method used for measuring the association between two
		    variables}
  \item{estimate }{a matrix of the correlation coefficient between variables}
  \item{pvalue }{a matrix of the p-value of the test}
  \item{nobs }{a matrix of the number of paired observation used in measuring the association between
	       pair of variables}
  \item{alternative }{a character string indicating the alternative hypothesis used}
}
\author{Alaine A. Gulles}
\keyword{univar}

