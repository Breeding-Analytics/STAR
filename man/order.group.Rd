\name{order.group}
\alias{order.group}
%- order.group.
\title{ Ordering the treatments according to the multiple comparison }
\description{
  This function allows us to compare the treatments averages or the adding of 
  their ranges with the minimal significant difference which can vary from one
  comparison to another one.
  This function is used by the HSD, LSD, Kruskal-Wallis, Friedman or Durbin 
  procedures.
}
\usage{
order.group(trt, means, N, MSerror, Tprob, std.err, parameter=1, snk=0, DFerror=NULL,alpha=NULL,sdtdif=NULL)
}

\arguments{
  \item{trt}{ Treatments }
  \item{means}{ Means of treatment }
  \item{N}{ Replications }
  \item{MSerror}{ Mean square error }
  \item{Tprob}{ minimum value for the comparison }
  \item{std.err}{ standard error }
  \item{parameter}{ Constante 1 (Sd), 0.5 (Sx) }
  \item{snk}{ Constante = 1 (Student Newman Keuls) }
  \item{DFerror}{ Degrees of freedom of the experimental error }
  \item{alpha}{ Level of risk for the test  } 
  \item{sdtdif}{ standar deviation of difference in BIB  } 
}
\value{
  \item{trt}{ Factor }
  \item{means}{ Numeric }
  \item{N}{ Numeric }
  \item{MSerror}{ Numeric }
  \item{Tprob}{ value between 0 and 1 }
  \item{std.err}{ Numeric }
  \item{parameter}{ Constant }
  \item{snk}{ Constant }
  \item{DFerror }{Numeric}
  \item{alpha }{Numeric}  
  \item{sdtdif}{Numeric}  
}

\author{ Felipe de Mendiburu }

\keyword{ manip }% at least one, from doc/KEYWORDS

