\name{tapply.stat}
\alias{tapply.stat}
%- tapply.stat.
\title{ Statistics of data grouped by factors }
\description{
  This process lies in finding statistics which consist of more than one variable,
  grouped or crossed by factors. The table must be organized by columns between 
  variables and factors.
}
\usage{
tapply.stat(y, x, stat = "mean")
}

\arguments{
  \item{y}{data.frame variables }
  \item{x}{data.frame factors  }  
  \item{stat}{ Method }
}

\value{
  \item{y}{ Numeric  }
  \item{x}{ Numeric }
  \item{stat}{method = "mean", ...  }
}
\author{ Felipe de Mendiburu }

\keyword{ univar }% at least one, from doc/KEYWORDS

