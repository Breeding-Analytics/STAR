FriedmanTest <- function(data, var, trmt, block) {
	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else { nameData <- paste(deparse(substitute(data))) }

	if (length(trmt) > 1) trmt <- trmt[1]
	if (length(block) > 1) block <- block[1]

	for (i in (1:length(var))) {
		if (any(table(data[,trmt], data[,block]) != 1)) {
			tempData <- aggregate(data[,var[i]], by = list(blk = data[,block], trt = data[,trmt]), FUN = mean, na.rm = TRUE)
			names(tempData) <- c(block, trmt, var[i])
			aggData <- TRUE
		} else {
			tempData <- data
			aggData <- TRUE
		}
		statResult <- data.frame(N = tapply(tempData[, var[i]],tempData[,trmt], length),
			     SUM.SCORE = tapply(tempData[, var[i]], tempData[,trmt], sum),
			     SUM.SCORE = tapply(tempData[, var[i]], tempData[,trmt], mean))
		statResult <- data.frame(rownames(statResult), statResult)
		colnames(statResult) <- c(trmt, "N", "Sum of Ranks", "Mean of Ranks")
		rownames(statResult) <- 1:nrow(statResult)

		cat("Rank Sums for Variable ", var[i], "\nClassified by ", trmt, "\n",sep = "")
		printDataFrame(statResult)
		cat("* Average scores are used for ties.\n\n")
		resultTest <- friedman.test(formula(paste(var, "~", trmt, "|", block)), tempData)

		if ((nchar(round(resultTest$statistic,0)) + 7) < nchar(resultTest$parameter)) { theWidth <- nchar(resultTest$parameter)
		} else { theWidth <- nchar(round(resultTest$statistic,0)) }
		cat(resultTest$method,"\n")
		cat(formatC(paste(rep("-", 27+theWidth), collapse = ""), width = 27+theWidth, format = "s"), "\n", sep = "")
		cat(formatC("Chi-Square", width = 20, format = "s", flag = "-"), formatC(resultTest$statistic, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
		cat(formatC("DF", width = 20, format = "s", flag = "-"), formatC(resultTest$parameter, width = theWidth + 7, format = "d"), "\n", sep = "")
		cat(formatC("Pr > Chi-Square", width = 20, format = "s", flag = "-"), formatC(resultTest$p.value, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
		cat(formatC(paste(rep("-", 27+theWidth), collapse = ""), width = 27+theWidth, format = "s"), "\n\n\n", sep = "")
	}
}