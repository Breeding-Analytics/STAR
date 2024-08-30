\name{LSD.test}
\alias{LSD.test}
%- LSD.
\title{ Multiple comparisons, "Least significant difference" and Adjust P-values }
\description{
  Multiple comparisons of treatments by means of LSD and a grouping of treatments. 
  The level by alpha default is 0.05. Returns p-values adjusted using one of several methods
}
\usage{
LSD.test(y, trt, DFerror, MSerror, alpha = 0.05, p.adj=c("none","holm", 
"hochberg", "bonferroni", "BH", "BY", "fdr"), group=TRUE, main = NULL, pwOrder = c("trmt","means"))
}

\arguments{
  \item{y}{ model(aov or lm) or answer of the experimental unit  }
  \item{trt}{ Constant( only y=model) or vector treatment applied to each experimental unit  }
  \item{DFerror}{ Degrees of freedom of the experimental error }
  \item{MSerror}{ Means square error of the experimental  }
  \item{alpha}{ Level of risk for the test  }
  \item{p.adj}{ Method for adjusting p values (see p.adjust)  }  
  \item{group}{ TRUE or FALSE }
  \item{main}{ title of the study  }
  \item{pwOrder}{means or trmt}
}  
\details{
  For equal or different repetition.
  p.adj = "holm", "hochberg", "bonferroni", "BH", "BY", "fdr".
  see p.adjust()
  p-adj ="none" is t-student.
  p-adj ="hommel" is not applied in this test.
}
\value{
 
  \item{y }{class (aov or lm) or vector numeric}
  \item{trt }{constant (only y=model) or vector alfanumeric}
  \item{DFerror }{Numeric}
  \item{MSerror }{Numeric}
  \item{alpha }{Numeric}
  \item{p.adj}{ text, see p.adjust  }  
  \item{group}{ Logic }
  \item{main }{Numeric}
}
\references{
     Steel, R.; Torri,J; Dickey, D.(1997)
     Principles and Procedures of Statistics
     A Biometrical Approach. pp178.
     }
\author{ Felipe de Mendiburu }

\seealso{\code{\link{HSD.test}}, \code{\link{SNK.test} }, \code{\link{duncan.test} } }

\keyword{ htest }% at least one, from doc/KEYWORDS

