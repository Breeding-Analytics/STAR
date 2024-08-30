# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# css
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

css <- function(x, na.rm = TRUE) UseMethod("css")

css.default <- function(x, na.rm = TRUE) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x) && !is.vector(x)) stop ("The argument should be a numeric vector.")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	css <- sum(x**2) - ((sum(x)**2)/length(x))
	return(css)
}

css.data.frame <- function(x, na.rm = TRUE) sapply(x, css)

