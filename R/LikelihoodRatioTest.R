#---------------------------------------------------------------------------------------
# R-CropStat Beta Version
# LikelihoodRatio.test Function: Functions for computing the likelihoodratio 
# Created by: Alaine A. Gulles 10.24.2010
# Modified by: Alaine A. Gulles 10.27.2011
# --------------------------------------------------------------------------------------

LikelihoodRatioTest <- function(x, y = NULL) UseMethod("LikelihoodRatioTest")

LikelihoodRatioTest.default <- function(x, y = NULL) {
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
	if (all(dim(x) == 2)) { result <- chisq.test(x, correct = TRUE) } else { result <- suppressWarnings(chisq.test(x, correct = FALSE)) }

	G <- 2 * sum(result$obs * log(result$obs/result$exp), na.rm = TRUE)
	pvalue <- 1 - pchisq(G, result[[2]][[1]])
	return(list(statistics = G, df = result[[2]][[1]], p.value = pvalue))
} ### end stmt -- LikelihoodRatioTest function
