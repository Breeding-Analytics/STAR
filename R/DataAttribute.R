# ----------------------------------------------------------------------------------------------------
# RCropStat Beta Version: Function
# ----------------------------------------------------------------------------------------------------
# DataAttribute
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.02.2012
# ----------------------------------------------------------------------------------------------------

DataAttribute <- function(data) UseMethod("DataAttribute")

DataAttribute.default <- function(data) {
	if(is.character(data)) { 
		nameData <- data
		if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
		tempData <- eval(parse(text = data))
	} else { 
		nameData <- paste(deparse(substitute(data)))	
		#if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
		tempData <- data
	}
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
	tempTable <- data.frame(names(tempData))
	withFactor <- FALSE
	for (i in (1:nrow(tempTable))) {
		if (is.factor(tempData[,i])) {
			withFactor <- TRUE
			if (is.ordered(tempData[,i])) tempTable[i,2] <- "ordered factor"
			else tempTable[i,2] <- "factor"
			tempTable[i,3] <- as.character(nlevels(tempData[,i]))
			if (nlevels(tempData[,i]) > 5) tempTable[i,4] <- paste(c(paste(levels(tempData[,i])[1:2], collapse = ", ", sep = ""),levels(tempData[,i])[nlevels(tempData[,i])]), collapse = ", ..., ", sep = "")
			else tempTable[i,4] <- paste(levels(tempData[,i]), collapse = ", ", sep = "")
		} else {
			if (typeof(tempData[,i]) == "double") { tempTable[i,2] <- "numeric" }
			else tempTable[i,2] <- typeof(tempData[,i])
			tempTable[i,3] <- ""
			tempTable[i,4] <- ""
		}
	}
	names(tempTable) <- c("VAR NAME", "TYPE", "NLEVELS", "LEVELS")
	if (!withFactor) { tempTable <- tempTable[,1:2] }
	return(tempTable)
}

