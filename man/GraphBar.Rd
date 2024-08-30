\name{GraphBar}
\alias{GraphBar}
\alias{GraphBar.default}
\title{GraphBar}
\description{
constructs bar graph
}
\usage{
GraphBar(data, outputPath, nVar, cVar, clustVars = NULL, mTitle = NULL, yAxisLab = NULL, 
	xAxisLab = NULL, minValue = NULL, maxValue = NULL, typeData = c("raw", "sumstat"), 
	descStat = c(NULL, "mean", "median", "sum", "freq"), barsHoriz = FALSE, barsClust = TRUE,
	byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdDev", "stdErr", "confInt"), 
	errMult = 1, confLev = NULL, axisLabelStyle = 1, barColor = NULL, 
	showLeg = FALSE, legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",
	"center"), legTitle = NULL, legCol = 1, boxed = TRUE,
	multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1,
	orientGraphs = c("left-right", "top-bottom"),
	showCatVarLevels = TRUE, barDensity, barLineAngle)
}
\arguments{
  \item{data}{a}
  \item{outputPath}{a}
  \item{nVar}{a}
  \item{cVar}{a}
  \item{clustVars}{a}
  \item{mTitle}{a}
  \item{yAxisLab}{a}
  \item{xAxisLab}{a}
  \item{minValue}{a}
  \item{maxValue}{a}
  \item{typeData}{a}
  \item{descStat}{a}
  \item{barsHoriz}{a}
  \item{barsClust}{a}
  \item{byVar}{a}
  \item{errBars}{a}
  \item{typeErrBar}{a}
  \item{errMult}{a}
  \item{confLev}{a}
  \item{axisLabelStyle}{a}
  \item{barColor}{a}
  \item{showLeg}{a}
  \item{legPos}{a}
  \item{legTitle}{a}
  \item{legCol}{a}
  \item{boxed}{a}
  \item{multGraphs}{a}
  \item{numRowsGraphs}{a}
  \item{numColsGraphs}{a}
  \item{orientGraphs}{a}
  \item{showCatVarLevels}{a}
  \item{barDensity}{a}
  \item{barLineAngle}{a}
}
\author{
Rose Imee Zhella A. Morantte
}
\keyword{bargraph}
