#---------------------------------------------------------------------------------------
# R-CropStat Beta Version
# ContingencyCoef Function: Functions for computing the contingency coefficient 
# Created by: Alaine A. Gulles 10.24.2010
# Modified by: Alaine A. Gulles 10.27.2011
# --------------------------------------------------------------------------------------

ContingencyCoef <- function(x, y = NULL) UseMethod("ContingencyCoef")

ContingencyCoef.default <- function(x, y = NULL) {
	DNAME <- deparse(substitute(x))
	if (is.data.frame(x)) x <- as.matrix(x)
	if (is.matrix(x)) { if (min(dim(x)) == 1) x <- as.vector(x)	}
	if (!is.matrix(x)) {
		if (is.null(y)) { stop("y must be a non-null vector") } else {
			if (length(x) != length(y)) { stop("x and y must have the same length") }
		}
		YDNAME <- deparse(substitute(y))
		ok <- complete.cases(x,y)
		x <- factor(x[ok])
		y <- factor(y[ok])
		if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) { stop("'x' and 'y' must have at least 2 levels") }
		x <- table(x,y)
		names(dimnames(x)) <- c(DNAME, YDNAME)
	}
	if (all(dim(x) == 2)) { 
		result <- chisq.test(x, correct = TRUE) 
		contingency.coef <- (prod(diag(result$obs)) - (result$obs[2,1]*result$obs[1,2]))/sqrt(prod(result$obs))
	} else { 
		result <- suppressWarnings(chisq.test(x, correct = FALSE)) 
		contingency.coef <- sqrt(result[[1]][[1]]/(sum(x) + result[[1]][[1]]))
	}
	return(contingency.coef)
} ### end stmt -- ContingencyCoef function
