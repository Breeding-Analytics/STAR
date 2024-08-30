# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# cv
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

cv <- function(x, na.rm = TRUE) UseMethod("cv")

cv.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	cv <- (sd(x)/mean(x))*100
	return(cv)
}

cv.data.frame <- function(x, na.rm = TRUE) sapply(x, cv)

coefVar <- function(x, na.rm = TRUE) cv(x, na.rm = TRUE)