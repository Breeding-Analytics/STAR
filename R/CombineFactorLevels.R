# --------------------------------------------------------------------------------------------
# CombineFactorLevels: Functions for combine levels of two or more factors into one variable
# Created by: Alaine A. Gulles 10.24.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 10.25.2011
# --------------------------------------------------------------------------------------------
# Arguments:
# data = a data frame or name of the data frame
# concatVar = names of the variables whose levels will be combine
# targetName = name of the new variable to be created
# --------------------------------------------------------------------------------------------

CombineFactorLevels <- function(data, concatVar, targetName = NULL) UseMethod("CombineFactorLevels")

CombineFactorLevels.default <- function(data, concatVar, targetName = NULL) {
	if (is.character(data)) { data <- eval(parse(text = data)) }
	if (is.null(targetName)) { tempName <- make.unique(c(names(data), "NewVar"), sep = "")
	} else { 
		if (!is.valid.name(targetName)) { targetName <- make.names(c(names(data), "NewVar"), unique = TRUE) }
		tempName <- make.unique(c(names(data), targetName), sep = "")
	}

	command <- paste("factor(paste(data[,'", paste(concatVar, collapse = "'],'_',data[,'", sep = ""),"'], sep = ''))", sep = "")
	data <- data.frame(data, eval(parse(text = command)))
	colnames(data)[ncol(data)] <- tempName[length(tempName)]
	remove(list = c("command", "tempName", "targetName"))
	return(data)
}
