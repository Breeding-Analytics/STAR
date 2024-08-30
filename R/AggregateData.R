# --------------------------------------------------------
# AggregateData : GUI for aggregating data sets
# Created by: Alaine A. Gulles 11.23.2010 
#		  for International Rice Research Institute
# Modified by: Alaine A. Gulles 03.28.2011 
# -------------------------------------------------------
# Arguments:
# data = name of a data frame or data frame
# var = names of the variable(s) to be aggregated
# grp = names of the grouping variable
# stat = statistics to be computed
# append = logical; is the aggregated values be append
#        with the original data
# -------------------------------------------------------

AggregateData <- function(data, var, grp, stat = c("mean"), append = FALSE) UseMethod("AggregateData")
	
AggregateData.default <- function(data, var, grp, stat = c("mean"), append = FALSE) {
	if (is.character(data)) { data <- eval(parse(text = data)) }
	# --- INPUT CHECKING:
     	#nameData <- paste(deparse(substitute(data)))    
	fxnLabel <- c("min", "max", "mean", "median","sum", "variance", "standard deviation")
	fxnEq <- c("min", "max", "mean", "median","sum", "var", "sd")
	
	capture.output(result <- DescriptiveStatistics(data, var, grp, statistics = fxnEq[match(stat, fxnLabel)]))
	aggrData <- reshape(result, idvar = grp,timevar = "Variable", direction = "wide")
	if (append) {
		tempData <- merge(data, aggrData, by = grp)
		tempData <- tempData[,unique(c(names(data), names(aggrData)))]
		remove(list = c("aggrData", "result", "fxnLabel", "fxnEq"))
		return(tempData)
	} else { 
		attr(aggrData, "reshapeWide") <- NULL
		remove(list = c("result", "fxnLabel", "fxnEq"))
		return(aggrData) 
	}
}

