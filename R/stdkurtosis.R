# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# stdkurtosis
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

stdkurtosis <- function(x, na.rm = TRUE) UseMethod("stdkurtosis")

stdkurtosis.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	stdkurtosis <- kurtosis(x)/sqrt(24/length(x))
	return(stdkurtosis)
}

stdkurtosis.data.frame <- function(x, na.rm = TRUE) sapply(x, stdkurtosis)

se.kurtosis <- function(x, na.rm = TRUE) stdkurtosis(x, na.rm = TRUE)