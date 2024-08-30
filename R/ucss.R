# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# ucss
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

ucss <- function(x, na.rm = TRUE) UseMethod("ucss")

ucss.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	ucss <- sum(x**2)
	return(ucss)
}

ucss.data.frame <- function(x, na.rm = TRUE) sapply(x, ucss)

