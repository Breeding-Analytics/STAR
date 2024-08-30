# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# stdskewness
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

stdskewness <- function(x, na.rm = TRUE) UseMethod("stdskewness")

stdskewness.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	stdskewness <- skew(x)/sqrt(6/length(x))
	return(stdskewness)
}

stdskewness.data.frame <- function(x, na.rm = TRUE) sapply(x, stdskewness)

se.skewness <- function(x, na.rm = TRUE) stdskewness(x, na.rm = TRUE)