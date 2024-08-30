# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# skewness
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

skewness <- function(x, na.rm = TRUE) UseMethod("skewness")

skewness.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	m3 <- sum((x - mean(x))^3)
	s3 <- sqrt(var(x))^3
	n <- length(x)
	skewness <- (n/((n-1)*(n-2)))*(m3/s3)
	return(skewness)
}

skewness.data.frame <- function(x, na.rm = TRUE) sapply(x, skewness)

skew <- function(x, na.rm = TRUE) skewness(x, na.rm = TRUE)