# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# kurtosis
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

kurtosis <- function(x, na.rm = TRUE) UseMethod("kurtosis")

kurtosis.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	m4 <- sum((x - mean(x))^4)
	s4 <- sqrt(var(x))^4
	n <- length(x)
	c4 <- (n * (n + 1))/((n - 1) * (n - 2) * (n - 3))
	kurtosis <- (c4 * (m4/s4)) - ((3 * (n - 1)**2)/((n - 2) * (n - 3)))
	return(kurtosis)
}

kurtosis.data.frame <- function(x, na.rm = TRUE) sapply(x, kurtosis)

coefKurtosis <- function(x, na.rm = TRUE) kurtosis(x, na.rm = TRUE)