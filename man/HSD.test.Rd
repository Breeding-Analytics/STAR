\name{HSD.test}
\alias{HSD.test}
%- HSD.
\title{ Multiple comparisons: Tukey  }
\description{
  It makes multiple comparisons of treatments by means of Tukey.
  The level by alpha default is 0.05.
}
\usage{
HSD.test(y, trt, DFerror, MSerror, alpha = 0.05, group=TRUE, main = NULL, pwOrder = c("trmt","means"))
}

\arguments{
  \item{y}{ model(aov or lm) or answer of the experimental unit  }
  \item{trt}{ Constant( only y=model) or vector treatment applied to each experimental unit  }
  \item{DFerror}{ Degree free }
  \item{MSerror}{ Mean Square Error }
  \item{alpha}{ Significant level }
  \item{group}{ TRUE or FALSE }
  \item{main}{ Title }
  \item{pwOrder}{means or trmt}
}
\details{
  It is necessary first makes a analysis of variance.
}
\value{

  \item{y }{class (aov or lm) or vector numeric}
  \item{trt }{constant (only y=model) or vector alfanumeric}
  \item{DFerror }{Numeric}
  \item{MSerror }{Numeric}
  \item{alpha }{Numeric}
  \item{group}{ Logic}  
  \item{main }{Text}

}
\references{ Principles and procedures of statistics a biometrical approach
             Steel & Torry & Dickey. Third Edition 1997 }
\author{ Felipe de Mendiburu author of package agricolae}

\seealso{\code{\link{LSD.test}}, \code{\link{SNK.test} },
\code{\link{duncan.test} } }
\keyword{ htest }

