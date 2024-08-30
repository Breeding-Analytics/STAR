\name{NormalityTest}
\alias{NormalityTest}
\alias{NormalityTest.default}
\title{Test for Normality}
\description{
  Performs test for normality using the Wilk-Shapiro, Shapiro-Francia,
  Cramer-von Mises, Anderson-Darling and Kolmogorov-Smirnov Test.
}
\usage{
  NormalityTest(data, var, grp = NULL, method = c("swilk"))
}
\arguments{
  \item{data}{a data frame or a character string}
  \item{var}{a non-empty character vector which indicates the name of the variable(s)}
  \item{grp}{NULL or a character string which indicates the grouping variable}
  \item{method}{a character string indicating which method to be used to test for normality}
}
\value{
  A data frame. 
}
\author{Alaine A. Gulles}
\keyword{normality}
\keyword{shapiro-Wilk}
\keyword{shapiro-Francia}
\keyword{Cramer-Von}
\keyword{Anderson}
\keyword{Kolmogorov-Smirnov}