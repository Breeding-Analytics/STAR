# -------------------------------------------------------------------------------
# R-CropStat Beta Version: 
# -------------------------------------------------------------------------------
# ConstructAOVTable: Functions for contructing the ANOVA table
# Created by: Alaine A. Gulles 05.18.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.18.2012
# -------------------------------------------------------------------------------
# Argument:
# aovTable is a class anova
# aovTable is a class "summary.aovlist" 
# aovTable is a class "summary.aov" and "listof"
# -------------------------------------------------------------------------------

ConstructAOVTable <- function(aovTable, computeTotal = TRUE) UseMethod("ConstructAOVTable")

ConstructAOVTable.default <- function(aovTable, computeTotal = TRUE) {
	temp <- aovTable
	rownames(temp)[nrow(temp)] <- "Error"
	if (computeTotal) {
		temp[nrow(temp)+1, 1:2] <- c(sum(temp[,1]), sum(temp[,2]))
		rownames(temp)[nrow(temp)] <- "Total"	
	}
	return(temp)

}

ConstructAOVTable.summary.aov <- function(aovTable, computeTotal = TRUE) {
	if (is.list(aovTable)) temp <- aovTable[[1]] else temp <- aovTable
	rownames(temp)[nrow(temp)] <- "Error"
	if (computeTotal) {
		temp[nrow(temp)+1, 1:2] <- c(sum(temp[,1]), sum(temp[,2]))
		rownames(temp)[nrow(temp)] <- "Total"	
	}
	return(temp)
}

ConstructAOVTable.summary.aovlist <- function(aovTable, computeTotal = TRUE) {
	
	temp <- NULL
	for (i in (1:length(aovTable))) {
		temp <- rbind(temp, aovTable[[i]][[1]])
		rownames(temp)[nrow(temp)] <- paste("Error(",letters[i],")", sep = "")
	}

	if (computeTotal) {
		temp[nrow(temp)+1, 1:2] <- c(sum(temp[,1]), sum(temp[,2]))
		rownames(temp)[nrow(temp)] <- "Total"	
	}
	return(temp)
} ### end stmt --- contructAOVTable function