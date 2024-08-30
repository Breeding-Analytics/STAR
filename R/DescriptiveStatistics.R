# ----------------------------------------------------------------------------------------
# RCropStat Beta Version: Function
# DescriptiveStatistics(data, var, grp = NULL, statToCompute)
# Created by: Alaine A. Gulles 10.01.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.04.2012 
# Note: Modified and combined version of the old DescriptiveStatistics and summaryStat
# ----------------------------------------------------------------------------------------

DescriptiveStatistics <- function(data, var, grp = NULL, statistics = c("nnmiss", "mean", "sd")) UseMethod("DescriptiveStatistics")

DescriptiveStatistics.default <- function(data, var, grp = NULL, statistics = c("nnmiss", "mean", "sd")) {
	if (is.character(data)) { 
		nameData <- data
		if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
		tempData <- eval(parse(text = data))
	} else {
		if (is.data.frame(data)) {
			nameData <- paste(deparse(substitute(data)))	
			tempData <- data
		} else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
	}
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
	if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
	if (!is.null(grp)) 		{ 
		if (any(is.na(match(grp, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
		tempGrp <- tempData[grp]
		for (i in (1:ncol(tempGrp))) { tempGrp[,i] <- factor(tempGrp[,i]) }
	} else { tempGrp <- rep(1, each = nrow(tempData)) }
	tempVar <- tempData[var]
	
	availableStatistics <- c("n", "nnmiss", "nmiss", "sum", "css", "ucss", "se.skew", "se.kurtosis", "range", "iqr", "var", "sd", "se.mean", 
		  "cv", "mean", "median", "mode", "min", "max", "q1", "q3", "skew", "kurtosis")
	if (is.null(statistics)) { statToCompute <- c("nnmiss", "mean", "sd") } else { statToCompute <- statistics }
	if (all(is.na(match(statToCompute, availableStatistics)))) { stop("All numerical descriptive measures enumerated is invalid. Expect one of the following: 'n', 'nnmiss', 'nmiss', 'sum', 'css', 'ucss', 'se.skew', 'se.kurtosis', 'range', 'iqr', 'var', 'sd', 'se.mean', 'cv', 'mean', 'median', 'mode', 'min', 'max', 'q1', 'q3', 'skew' or 'kurtosis'") }
	statToCompute <- availableStatistics[na.omit(match(statToCompute, availableStatistics))]

	summaryTable <- NULL
	outputLabel <- NULL
	
	# --- compute for number of rows in the data set
	if(!is.na(match("n", statToCompute))) { 
		a <- NULL
		for (i in (1:ncol(tempVar))) { a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, length))) }
		if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		summaryTable <- a
		colnames(summaryTable)[ncol(summaryTable)] <- "N_Obs"
		outputLabel <- c(outputLabel, "No. of Obs")
	}

	# --- compute for the number of non-missing observations
	if(!is.na(match("nnmiss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) {
			newData <- na.omit(cbind(tempVar[i], tempGrp))
			a <- rbind(a, as.data.frame.table(tapply(newData[,1], newData[2:ncol(newData)], length)))
		}
		if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "N_NonMissObs"
		outputLabel <- c(outputLabel, "No. of Non-Missing Obs.")
	}

	# --- compute for the number of missing observations
	if(!is.na(match("nmiss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) {
			newData <- cbind(tempVar[i], tempGrp)
			newData <- newData[!complete.cases(newData),]
			if (!is.null(grp)) a <- rbind(a, as.data.frame.table(tapply(newData[,1], newData[2:ncol(newData)], length)))
			else a <- rbind(a, data.frame(tempGrp = 1, Freq = nrow(newData)))
		}
		a$Freq <- replace(a$Freq, is.na(a[,"Freq"]), 0)
		#if (any(is.na(a$Freq))) { a$Freq <- replace(a$Freq, attr(na.omit(a$Freq), "na.action"), 0) }
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "N_MissObs"
		outputLabel <- c(outputLabel, "No. of Missing Obs.")
	}
	
	# --- compute the minimun observation
	if(!is.na(match("min", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Min"
		outputLabel <- c(outputLabel, "Minimum")
	}

	# --- compute the maximum observation
	if(!is.na(match("max", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, max, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Max"
		outputLabel <- c(outputLabel,"Maximum")
	}

	# --- compute the sum of the variable
	if(!is.na(match("sum", statToCompute))) {
		a <- NULL		
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, sum, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Sum"
		outputLabel <- c(outputLabel,"Sum")
	}

	# --- compute the mean of the variable
	if(!is.na(match("mean", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, mean, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Mean"
		outputLabel <- c(outputLabel,"Mean")
	}

	# --- compute the median of the variable
	if(!is.na(match("median", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, median, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Median"
		outputLabel <- c(outputLabel,"Median")
	}
	# --- determine the modal value of the variable
	if(!is.na(match("mode", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, data.frame(Var = i, Freq = paste(tapply(tempVar[[i]], tempGrp, modalValue, na.rm = TRUE)[[1]], collapse = ", ", sep = "")))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Mode"
		outputLabel <- c(outputLabel,"Mode")
	}
	if(!is.na(match("q1", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q1"
		outputLabel <- c(outputLabel,"1st Quartile")
	}
	if(!is.na(match("q3", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, quantile, probs = 0.75, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q3"
		outputLabel <- c(outputLabel,"3rd Quartile")
	}
	if(!is.na(match("range", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, max, na.rm = TRUE) - tapply(tempVar[[i]], tempGrp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Range"
		outputLabel <- c(outputLabel,"Range")
	}
	if(!is.na(match("iqr", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, quantile,probs = 0.75, na.rm = TRUE) - tapply(tempVar[[i]], tempGrp, quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "IQR"
		outputLabel <- c(outputLabel,"Inter Quartile Range")
	}
	if(!is.na(match("var", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, sd, na.rm = TRUE)**2))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Variance"
		outputLabel <- c(outputLabel,"Variance")
	}
	if(!is.na(match("sd", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, sd, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "StdDev"
		outputLabel <- c(outputLabel,"Standard Deviation")
	}
	if(!is.na(match("se.mean", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stdmean, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Mean"
		outputLabel <- c(outputLabel,"Std. Error of the Mean")
	}
	if(!is.na(match("cv", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, cv, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "CV"
		outputLabel <- c(outputLabel,"Coefficient of Variation")
	}
	if(!is.na(match("css", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, css, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "CSS"
		outputLabel <- c(outputLabel,"Corrected Sum of Squares")
	}
	if(!is.na(match("ucss", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, ucss, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "UCSS"
		outputLabel <- c(outputLabel,"Uncorrected Sum of Squares")
	}
	if(!is.na(match("skew", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, skewness, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Skewness"
		outputLabel <- c(outputLabel,"Skewness")
	}
	if(!is.na(match("se.skew", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stdskewness, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Skew"
		outputLabel <- c(outputLabel,"Std. Error of Skewness")
	}
	if(!is.na(match("kurtosis", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, kurtosis, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Kurtosis"
		outputLabel <- c(outputLabel,"Kurtosis")
	}
	if(!is.na(match("se.kurtosis", statToCompute))) {
		a <- NULL
		for (i in (1:ncol(tempVar))) a <- rbind(a, as.data.frame.table(tapply(tempVar[[i]], tempGrp, stdkurtosis, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE_Kurtosis"
		outputLabel <- c(outputLabel,"Std. Error of Kurtosis")
	}
	
	if (is.null(ncol(tempGrp))) { summaryTable[,1] <- names(tempVar)
	} else {
		variable <- c(rep(names(tempVar), each = nrow(as.data.frame.table(table(tempGrp)))))
		summaryTable <- data.frame(variable,summaryTable)
	}
	colnames(summaryTable)[1] <- "Variable"

	options(width = 5000)
	cat("Descriptive Statistics\n")
	printDataFrame(summaryTable)
	return(invisible(summaryTable))
}
