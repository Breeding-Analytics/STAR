# -------------------------------------------------------------------------------
# STAR
# -------------------------------------------------------------------------------
# PredictedValues
# Created by: Alaine A. Gulles for International Rice Research Institute
# -------------------------------------------------------------------------------


PredictedValues <- function(object, errorTerm = NULL) UseMethod("PredictedValues")

PredictedValues.default <- function(object, errorTerm = NULL) {
	return(fitted(object))
}

PredictedValues.aovlist <- function(object, errorTerm = NULL) {
	if (!inherits(object, what = "aovlist")) stop("Not an aovlist object.")
	aov.proj <- proj(object)
	if (is.null(errorTerm)) { numGrp <- length(object)
	} else { numGrp <- which(names(proj(object)) == errorTerm)}
	resultTable <- proj(object)[[1]][,1]
	for (i in (2:numGrp)) {
		nterms <- ncol(proj(object)[[i]])
		if (dimnames(proj(object)[[i]])[[2]][nterms] == "Residuals") nterms <- nterms - 1
		if (nterms > 0) {
			if (nterms == 1) { resultTable <- resultTable + proj(object)[[i]][,1]
			} else { resultTable <- resultTable + rowSums(proj(object)[[i]][,1:nterms]) }
		}
	}
	return(resultTable)
}

