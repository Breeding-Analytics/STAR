# ----------------------------------------------------------------------------------
# ToLong: Function for restructuring data from wide to long format
# Created By: Alaine A. Gulles 11.26.2010 for International Rice Research Institute
# Modified By: Alaine A. Gulles 07.09.2012
# ----------------------------------------------------------------------------------

ToLong <- function(data, varying, timevar, vnames, idvar, label = NULL) UseMethod("ToLong")

ToLong.default <- function(data, varying, timevar, vnames, idvar, label = NULL) {

	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else {
		if (!is.data.frame(data)) { stop("The argument 'data' should be a data frame.") }
		nameData <- paste(deparse(substitute(data)))
	}
	vnames <- unique(vnames)

	if (is.list(varying)) {
		tempVarying <- unlist(varying)
		lengthVar <- length(varying[[1]])
		tempVar <- paste("c('", paste(varying[[1]], collapse = "', '", sep = ""),"')", sep = "")
		if (length(varying) > 1) {
			for (i in (2:length(varying))) {
				tempVar <- paste(tempVar, paste("c('", paste(varying[[i]], collapse = "', '", sep = ""),"')", sep = ""), sep = ", ")
				lengthVar <- c(lengthVar, length(varying[[i]]))
			}
		}
		if (length(varying) != length(vnames)) {
			if (length(varying) > length(vnames)) { vnames <- make.unique(c(vnames, paste(rep("newVar", (length(varying)-length(vnames))), 1:(length(varying)-length(vnames)),sep = "")), sep = "")
			} else {  vnames <- vnames[1:length(varying)] } 
		}
	} else {
		tempVarying <- varying
		lengthVar <- length(tempVarying)
		tempVar <- paste("c('", paste(varying, collapse = "','", sep = ""),"')", sep = "")
	}

	if (!is.null(label)) { if (length(label) != max(lengthVar)) { label <- c(1:max(lengthVar)) }
	} else { label <- c(1:max(lengthVar)) }
     
     # determine if all elements of a component is of the same data structure
     # if not, convert all elements to a factor
     for (i in (1:length(varying))) {
          if (any(DataAttribute(data[,varying[[i]]])[2] == "factor")) {
               if (!all(DataAttribute(data[,varying[[i]]])[2] == "factor")) {
                    varTemp <- as.character(DataAttribute(data[,varying[[i]]])[which(DataAttribute(data[,varying[[i]]])[2] != "factor"),1])
                    for (j in (1:length(varTemp))) {
                         data[,varTemp[j]] <- factor(data[,varTemp[j]])
                    } 
               }
          }
     }
     
     #if (length(unique(data[,idvar])) == length(data[,idvar])) {
	#     #command <- paste("reshape(data = ",nameData,", v.names = c('",paste(vnames, collapse = "', '", sep = ""),"')", sep = "")
	#     command <- paste("reshape(data = data, v.names = c('",paste(vnames, collapse = "', '", sep = ""),"')", sep = "")
	#     command <- paste(command, ", varying = list(", tempVar,")", sep = "")
	#     command <- paste(command, ", idvar = c('", paste(idvar,collapse = "','", sep = ""),"')", sep = "")
	#     command <- paste(command, ", timevar = c('", paste(timevar, collapse = "', '", sep = ""),"')", sep = "")
	#     
	#     # determine if there are variable to be drop
	#     if (length(c(tempVarying, idvar)) == length(names(data))) {
	#          command <- paste(command, ", drop = NULL", sep = "")
	#     } else {
	#          index <- match(c(tempVarying, idvar), names(data))
	#          varToDrop <- names(data)[-I(index)]
	#          command <- paste(command, ", drop = c('", paste(varToDrop, collapse = "', '", sep = ""),"')", sep = "")
	#          remove(list = c("index", "varToDrop"))
	#     }
	#     
	#     command <- paste(command, ", times = c('", paste(label, collapse = "', '", sep = ""),"'), direction = 'long')", sep = "")
	#     newData <- eval(parse(text = command))
	#     
	#} else {
          newData <- NULL
          for (i in (1:length(varying[[1]]))) {
               commandSub <- NULL
               for (j in (1:length(varying))) {
                    if (is.null(commandSub)) {
                         commandSub <- paste("'",varying[[j]][i], "'", sep = "")
                    } else {
                         commandSub <- paste(commandSub, paste(", '",varying[[j]][i], "'", sep = ""), sep = "")     
                    }
               }
               command <- paste("data[,c('", paste(idvar, collapse = "', '"),"', ", commandSub,")]", sep = "")
               tempData <- eval(parse(text = command))
               names(tempData)[(ncol(tempData) - length(varying) + 1):ncol(tempData)] <- vnames
               newData <- rbind(newData, data.frame(label[i], tempData))
          }
          newData <- newData[,c(idvar, names(newData)[1], vnames)]
          names(newData)[length(idvar)+1] <- timevar
	#}
     

	rownames(newData) <- 1:nrow(newData)
	for (i in (1:length(timevar))) { newData[,timevar[i]] <- factor(newData[,timevar[i]]) }
	attr(newData, "reshapeLong") <- NULL
	remove(list = c("nameData", "tempVar", "lengthVar", "tempVarying", "command"))
	return(newData)
}