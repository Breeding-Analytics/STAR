KruskalWallisTest <- function(data, var, grp) UseMethod("KruskalWallisTest")

KruskalWallisTest.default <- function(data, var, grp) {

	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else { nameData <- paste(deparse(substitute(data))) }

	resultTable <- NULL
	for (i in (1:length(grp))) {
		for (j in (1:length(var))) {
			resultStat <- data.frame(N = tapply(data[,var[j]], data[,grp[i]], length),
							 SUM.SCORE = tapply(rank(data[,var[j]]), data[,grp[i]], sum),
							 MEAN.SCORE = tapply(rank(data[,var[j]]), data[,grp[i]], mean))
			resultStat <- data.frame(rownames(resultStat), resultStat)
			colnames(resultStat) <- c(grp[i], "N", "Sum of Ranks", "Mean of Ranks")
			rownames(resultStat) <- 1:nrow(resultStat)

			cat("Rank Sums for Variable ", var[j],"\nClassified by Variable ",grp[i],"\n", sep = "")
			printDataFrame(resultStat)
			cat("* Average scores are used for ties.\n\n\n")
			resultTest <- kruskal.test(data[,var[j]], data[,grp[i]])
			if ((nchar(round(resultTest$statistic,0)) + 7) < nchar(resultTest$parameter)) { theWidth <- nchar(resultTest$parameter)
			} else { theWidth <- nchar(round(resultTest$statistic,0)) }
			cat(resultTest$method,"\n")
			cat(formatC(paste(rep("-", 27+theWidth), collapse = ""), width = 27+theWidth, format = "s"), "\n", sep = "")
			cat(formatC("Chi-Square", width = 20, format = "s", flag = "-"), formatC(resultTest$statistic, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
			cat(formatC("DF", width = 20, format = "s", flag = "-"), formatC(resultTest$parameter, width = theWidth + 7, format = "d"), "\n", sep = "")
			cat(formatC("Pr > Chi-Square", width = 20, format = "s", flag = "-"), formatC(resultTest$p.value, width = theWidth + 7, format = "f", digits = 4), "\n", sep = "")
			cat(formatC(paste(rep("-", 27+theWidth), collapse = ""), width = 27+theWidth, format = "s"), "\n\n\n", sep = "")

			resultTable <- rbind(resultTable, data.frame(Group = grp[i], Variable = var[j], Statistics = resultTest$statistic, DF = resultTest$parameter, pvalue = resultTest$p.value))
		}
	}
	rownames(resultTable) <- 1:nrow(resultTable)
	return(invisible(resultTable))
}