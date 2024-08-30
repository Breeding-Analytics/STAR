# --------------------------------------------------------------------------------------
# R-CropStat Beta Version: Function
# --------------------------------------------------------------------------------------
# foldedFTest: Perform F-test to compare variances of two samples
#              from normal populations. Uses the folded form of the
#		   two-tailed F statistic, which test the hypothesis that 
#		   the variances of two populationsare equal. 
# Created by: Alaine A. Gulles 02.24.2011 for International Rice Research Institute
# Modified by : Alaine A. Gulles 02.24.2011
# --------------------------------------------------------------------------------------
# ARGUMENT: x   = a numeric vector of data values
#           grp = a grouping vector, can be numeric or character 
#
# --------------------------------------------------------------------------------------

FoldedFTest <- function(x, grp) {

	result <- rbind(tapply(x, grp, var), tapply(x, grp, length))
	result <- result[,order(result[1,])]

	fc = result[1,ncol(result)][[1]]/result[1,1][[1]]
	pval <- 2*pf(fc, result[2,ncol(result)][[1]] - 1, result[2,1][[1]]- 1, lower.tail = FALSE) 
	df <- c(NUMDF = result[2,ncol(result)][[1]]- 1, DENDF = result[2,1][[1]]- 1)
	return(list(statistic = fc[[1]], 
      	    df = c(num.df = result[2,ncol(result)][[1]]- 1, denom.df = result[2,1][[1]]- 1), 
	          pvalue = pval[[1]],
		    method = "Folded F"))

}### end stmt --- Foldedf.test