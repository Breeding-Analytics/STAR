# -----------------------------------------------------------------------------------------------------------------
# SortCases: Rearrange the observations based on the value of one or more sorting variables.
# Created by: Alaine A. Gulles 10.15.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 08.31.2012
# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------
# Arguments: 
# data = name of the data frame or a data frame
# var = a non-empty character string or vector which indicate the names of the variable(s) that will be 
#	the basis for sorting the data frame 
# sortBy = a character string or vector which indicates the order by which the variable(s) will be sorted
#	= possible values are "ascending" and/or "descending"
# -----------------------------------------------------------------------------------------------------------------

SortCases <- function(data, var, sortBy = "ascending") UseMethod("SortCases")

SortCases.default <- function(data, var, sortBy = "ascending") {
	if (is.character(data)) { data <- eval(parse(text = data))	}
     	if (length(sortBy) != length(var)) {
         	if (length(sortBy) > length(var)) { sortBy <- sortBy[1:length(var)] 
          	} else {
            	if (length(sortBy) > 1) { sortBy <- c(sortBy, rep("ascending", length(var) - length(sortBy)))   
               	} else { sortBy <- rep(sortBy, length(var)) }
          	}
     	}

     	for (i in (1:length(var))) {
          	if (rev(sortBy)[i] ==  "ascending") {
               	data <- data[order(data[,rev(var)[i]], decreasing = FALSE),]
          	} else { data <- data[order(data[,rev(var)[i]], decreasing = TRUE),] }
     	}
     	rownames(data) <- 1:nrow(data)
     	return(data)     
}