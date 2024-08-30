# -------------------------------------------------------------------------------
# STAR
# -------------------------------------------------------------------------------
# ClassInformation
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles and Rose Imee Zhella A. Morantte 04.16.2013
# Note: Include multiple response variable
# -------------------------------------------------------------------------------

ClassInformation <- function(data, respvar = NULL) UseMethod("ClassInformation")

ClassInformation.default <- function(data, respvar = NULL) {
	result <- DataAttribute(data)
	result <- result[result[,"TYPE"] == "factor",c(1,3,4)]
	result[,2] <- as.numeric(result[,2])
	result[,2] <- factor(result[,2])
	colnames(result) <- c("FACTOR", "NO. OF LEVELS", "LEVELS")
	#cat("Class Level Information\n")
	cat("Summary Information\n")
	printDataFrame(result)
	#cat("\n")
	if (!is.null(respvar) && length(respvar) == 1) {
		if (length(na.omit(data[,respvar])) == nrow(data)) {
			cat("Number of Observations Read and Used:",nrow(data),"\n")
		} else {
			cat("Number of Observations Read:",nrow(data),"\n")
			cat("Number of Observations Used:",length(na.omit(data[,respvar])),"\n\n")
		}
	}
	if (!is.null(respvar) && length(respvar) > 1) {
		numObsNM = nrow((na.omit(data[,respvar])))
		numObs = nrow((data[,respvar]))

		if (numObs == numObsNM) { cat("Number of Observations:", numObs,"\n")
		} else {
			cat("Number of Observations:", numObs,"\n")
			cat("Number of Observations Used:", numObsNM,"\n\n")
		}

	}
}

