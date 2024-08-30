# ----------------------------------------------------------------------------------
# Restructure: function for restructuring data from long to wide format
# Created By: Alaine A. Gulles 01.26.2011 for International Rice Research Institute
# Modified By: Alaine A. Gulles 07.10.2012
# ----------------------------------------------------------------------------------

ToWide <- function(data, vnames, timevar, idvar) UseMethod("ToWide")

ToWide.default <- function(data, vnames, timevar, idvar) {
	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else { 
		if (!is.data.frame(data)) { stop("data should be a data frame.") }
		nameData <- paste(deparse(substitute(data)))
	}

	# determine if there are variables to be drop
	if (length(c(vnames, idvar, timevar)) == length(names(data))) {  varToDrop <- "NULL"
	} else {
		index <- match(c(vnames, idvar, timevar), names(data))
		varToDrop <- names(data)[-I(index)]	
		varToDrop <- paste("c('", paste(varToDrop, collapse = "', '"),"')", sep = "")
	}

	if (length(timevar) == 1) {
		command <- paste("reshape(data = ",nameData,", v.names = c('",paste(vnames, collapse = "', '"),"'),", sep = "") 
		command <- paste(command, " idvar = c('", paste(idvar, collapse = "', '"),"'),", sep = "")
		command <- paste(command, " timevar = c('", paste(timevar, collapse = "', '"),"'),", sep = "")
		command <- paste(command, " drop = ", varToDrop,", direction = 'wide')", sep = "")
		newData <- eval(parse(text = command))
		#if (nrow(newData) == 0) {  stop("The restructure data is empty.") }
	} else {
		timevartemp1 <- timevar[1]
		timevartemp2 <- timevar[2:length(timevar)]
		idvartemp1 <- c(idvar, timevartemp2)
		newData <- NULL
		tempData <- data
		while(!is.na(timevartemp2)) {
			if (!is.null(newData)) {
				timevartemp1 <- timevartemp2[1]
				timevartemp2 <- timevartemp2[2:length(timevar)]
				if (is.na(timevartemp2)) { idvartemp1 <- c(idvar) 
				} else { idvartemp1 <- c(idvar, timevartemp2) }
				index <- match(c(timevartemp1, idvartemp1),names(newData))
				vnames <- names(newData[-I(index)])
				tempData <- newData
				varToDrop <- "NULL"
			}
			command <- paste("reshape(data = tempData, v.names = c('",paste(vnames, collapse = "', '"),"'),", sep = "") 
			command <- paste(command, " idvar = c('", paste(idvartemp1, collapse = "', '"),"'),", sep = "")
			command <- paste(command, " timevar = c('", paste(timevartemp1, collapse = "', '"),"'),", sep = "")
			command <- paste(command, " drop = ", varToDrop,", direction = 'wide')", sep = "")
			newData <- eval(parse(text = command))
		}
		remove(list = c("timevartemp1", "timevartemp2", "idvartemp1", "tempData", "index"))
	}
	remove(list = c("varToDrop", "command"))
	rownames(newData) <- 1:nrow(newData)
	attr(newData, "reshapeWide") <- NULL
	return(newData)
}