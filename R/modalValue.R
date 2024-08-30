# ------------------------------------------------------------------------
# RCropStat Beta Version: Function for Utilities Statistics
# ------------------------------------------------------------------------
# modalValue
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.03.2012
# ------------------------------------------------------------------------

modalValue <- function(x, na.rm = TRUE) UseMethod("modalValue")

modalValue.default <- function(x, na.rm = TRUE) {
	if (!is.vector(x)) { stop ("The argument should be a numeric vector.") }
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	f <- table(x)
	f <- sort(f)
	if (length(f) == 1) mode <- x[1]
	else {
		if (f[[1]] == f[[length(f)]]) mode <- ""
		else {
			i <- length(f) - 1
			modalValue <- as.character(rownames(f)[length(f)])
			while(f[[i]] == f[[length(f)]]) {
				modalValue <- c(modalValue,rownames(f)[i])
				i <- i-1
			}
			modalValue <- sort(modalValue)	
		}
	}
	return(modalValue)
}


modalValue.data.frame <- function(x, na.rm = TRUE) sapply(x, modalValue)

