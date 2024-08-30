#----------------------------------------------------------------
# getVarInfo                        
# returns a data frame containing names and corresponding types 
# of variables in a data set
#
# ARGUMENTS:
# fileName - name of the input file
# fileFormat - 1 if input file is of .csv format;
#		   2 if it is a .txt file; or
#		   3 if it is a .rda file
# separator - for .txt files; separator used in the file
#----------------------------------------------------------------

getVarInfo <- function(fileName, fileFormat, separator = NULL) {
	
	#Read file
	#if .csv
	if (fileFormat == 1) {
		data <- read.csv(fileName,header=TRUE, na.strings = c("NA","."), blank.lines.skip=TRUE, sep = ",")
	}
	#if .txt 
 	else if (fileFormat == 2) {
		if (is.na(separator)) separator = "\t"
		data <- read.table(fileName,header=TRUE, na.strings = c("NA","."),blank.lines.skip = TRUE, sep = separator)

	#if .rda
	} else {
		data <- eval(parse(text = load(fileName)))
	}

	#save names of variables into a data frame
	varInfo = as.data.frame(names(data))

	for (i in 1:ncol(data)) {
		if (is.numeric(data[,i]))
			varInfo[i,2] <- "Numeric"
		else varInfo[i,2] <- "Factor"
	}

	names(varInfo) = c("Variable", "Type")
	
	return(varInfo)
}