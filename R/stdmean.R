# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# stdmean
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

stdmean <- function(x, na.rm = TRUE) UseMethod("stdmean")

stdmean.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	stdmean <- sqrt(var(x)/length(x))
	return(stdmean)
}

stdmean.data.frame <- function(x, na.rm = TRUE) sapply(x, stdmean)

se.mean <- function(x, na.rm = TRUE) stdmean(x, na.rm = TRUE)