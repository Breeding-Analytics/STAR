# -------------------------------------------------------------------------------------
# R-SCRIPT FOR STAR:
# GoodnessOfFitTest: Perform Chi-Square Goodness of Fit test
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 
# -------------------------------------------------------------------------------------

GoodnessOfFitTest <- function(data, var, colFreq = NULL, expected = NULL) UseMethod("GoodnessOfFitTest")
     
GoodnessOfFitTest.default <- function(data, var, colFreq = NULL, expected = NULL) {

	if (is.character(data)) { data <- eval(parse(text = data)) } 
	if (is.null(colFreq)) { tempObs <- table(data[,var]) 
	} else {
		if (length(data[,var]) != length(data[,colFreq])) { stop("The argument 'var' and 'colFreq' should be of the same length.")}
		tempObs <- as.table(data[,c(colFreq)])
		names(tempObs) <- data[,var]
	}
	if (is.null(expected)) {
		tempExp <- tempObs
		for (i in (1:length(tempExp))) { tempExp[i] <- sum(tempObs)/length(tempObs) }
	} else {
		tempExp <- expected
		tempExp <- tempExp[order(names(tempExp))]
	}

	result <- chisq.test(tempObs, p = tempExp, rescale.p = TRUE, correct = FALSE)
    	tempTable <- cbind(result$observed, result$expected)
    	tempTable <- data.frame(row.names(tempTable), tempTable)
    	colnames(tempTable) <- c(var, "Obs Freq", "Exp Freq")
    	rownames(tempTable) <- 1:nrow(tempTable)
    	cat("Frequency Table:\n")
    	printDataFrame(tempTable)
   	cat("\n\n")
    	if ((nchar(round(result$statistic, 0)) + 7) < nchar(result$parameter)) {
      	theWidth <- nchar(result$parameter)
    	} else { theWidth <- nchar(round(result$statistic, 0)) }
    	cat("Chi-Square Goodness of Fit Test\n")
    	cat(formatC(paste(rep("-", 27 + theWidth), collapse = ""), width = 27 + theWidth, format = "s"), "\n", sep = "")
    	cat(formatC("Chi-Square", width = 20, format = "s", flag = "-"), formatC(result$statistic, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
    	cat(formatC("DF", width = 20, format = "s", flag = "-"), formatC(result$parameter, width = theWidth + 7, format = "d"), "\n", sep = "")
    	cat(formatC("Pr > Chi-Square", width = 20, format = "s", flag = "-"), formatC(result$p.value, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
    	cat(formatC(paste(rep("-", 27 + theWidth), collapse = ""), width = 27 + theWidth, format = "s"), "\n\n\n", sep = "")
	return(invisible(list(FreqTable = tempTable, statistics = result$statistic, df = result$parameter, pvalue = result$p.value)))	
}