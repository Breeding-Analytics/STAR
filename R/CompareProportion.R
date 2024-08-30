# -------------------------------------------------------------------------------------
# R-SCRIPT FOR STAR:
# CompareProportion: Compare proportion for one and two population.
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 
# -------------------------------------------------------------------------------------

CompareProportion <- function(data, varX, varY = NULL, 
	grp = NULL, testVal = NULL, procedure = c("one", "paired", "independent"),
	alternative = c("two.sided", "greater", "less"),
	CI = FALSE, confLevel = 0.95) UseMethod("CompareProportion")

CompareProportion.default <- function(data, varX, varY = NULL, 
	grp = NULL, testVal = NULL, procedure = c("one", "paired", "independent"),
	alternative = c("two.sided", "greater", "less"),
	CI = FALSE, confLevel = 0.95) {

	if (is.character(data)) { 
		nameData <- data
		data <- eval(parse(text = data))
	} else {
		if (is.data.frame(data)) { nameData <- paste(deparse(substitute(data)))	
		} else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
	}

	alternative <- match.arg(alternative)
	procedure <- match.arg(procedure)
	
	if (procedure == "one") {
		estResult <- NULL
		summaryResult <- NULL
		#if (length(successLevel) != length(varX)) { stop("varX and successLevel should be of the same length.")}
		if (is.null(testVal)) {
			testVal <- 0.5
		} 
		if (length(testVal) != length(varX)) {
			if (length(testVal) < length(varX)) { testVal <- c(testVal, rep(0.5,(length(varX) - length(testVal))))
			} else { testVal <- testVal[1:length(varX)] }
		}


		for (i in (1:length(varX))) {
			tempTable <- table(data[, varX[i]])
			if (length(tempTable) != 2) { 
                    cat(paste("ERROR: The variable '",varX[i],"' should be a dichotomous variable.\n\n", sep = "")) 
                    next
			} else { 
				if (any(is.na(match(names(tempTable), c("0","1"))))) {
					cat(paste("ERROR: The variable '",varX[i],"' should be a dichotomous variable with value 1 for success and 0 for failure.\n\n", sep = ""))
                         next
				} else { names(tempTable) <- c("Failure", "Success") }
			}
			successCnt <- tempTable["Success"]
			totalCnt <- sum(tempTable)

			chiResult <- prop.test(successCnt, totalCnt, p = testVal[i], conf.level = confLevel, correct = FALSE)
			exactResult <- binom.test(successCnt, totalCnt, p = testVal[i], conf.level = confLevel)
			approxZ <- ((successCnt/totalCnt) - testVal[i])/sqrt((testVal[i]*(1-testVal[i]))/totalCnt)
			if (alternative == "less") {	approxZ.pval <- pnorm(approxZ)
			} else {
				if (alternative == "greater") { approxZ.pval <- pnorm(approxZ, lower.tail = FALSE)
				} else { approxZ.pval <- 2*pnorm(-abs(approxZ)) }
			}
		
			if (i == 1) { 
				if (CI) {
					approxZWidth <-  qnorm(1-((1-confLevel)/2))*sqrt((testVal[i]*(1-testVal[i]))/totalCnt)
					estResult <- data.frame(Variable = varX[i], SuccessCnt = successCnt, 
									Total = totalCnt, Estimate = chiResult$estimate,
									LL1 = chiResult$conf.int[1], UL1 = chiResult$conf.int[2], 
									LL2 = exactResult$conf.int[1], UL2 = exactResult$conf.int[2],
									LL3 = (successCnt/totalCnt) - approxZWidth, UL3 = (successCnt/totalCnt) + approxZWidth)
				} else { estResult <- data.frame(Variable = varX[i], SuccessCnt = successCnt, Total = totalCnt, Estimate = chiResult$estimate)	}
				summaryResult <- data.frame(Variable = varX[i],Procedure = "Chi-Square", TestValue = testVal[i], Statistics = chiResult$statistic, Prob = chiResult$p.value)
			} else {
				if (CI) { 
					approxZWidth <-  qnorm(1-((1-confLevel)/2))*sqrt((testVal[i]*(1-testVal[i]))/totalCnt)
					estResult <- rbind(estResult,
							 data.frame(Variable = varX[i], SuccessCnt = successCnt, 
									Total = totalCnt, Estimate = chiResult$estimate,
									LL1 = chiResult$conf.int[1], UL1 = chiResult$conf.int[2], 
									LL2 = exactResult$conf.int[1], UL2 = exactResult$conf.int[2],
									LL3 = (successCnt/totalCnt) - approxZWidth, UL3 = (successCnt/totalCnt) + approxZWidth))
				} else { estResult <- rbind(estResult, data.frame(Variable = varX[i], SuccessCnt = successCnt, Total = totalCnt, Estimate = chiResult$estimate))	}
				summaryResult <- rbind(summaryResult,data.frame(Variable = varX[i],Procedure = "Chi-Square", TestValue = testVal[i], Statistics = chiResult$statistic, Prob = chiResult$p.value))
			} 
			summaryResult <- rbind(summaryResult, data.frame(Variable = varX[i], Procedure = "Exact Binomial", TestValue = testVal[i], Statistics = exactResult$statistic, Prob = exactResult$p.value))
			summaryResult <- rbind(summaryResult, data.frame(Variable = varX[i], Procedure = "Approx Z", TestValue = testVal[i], Statistics = approxZ, Prob = approxZ.pval))
		}
		rownames(summaryResult) <- 1:nrow(summaryResult)
		if (CI) { colnames(estResult) <- c("Variable", "Num Success", "Num Trial", "Proportion", 
						 "ChiSq LL*", "ChiSq UL*", "Binomial LL**", "Binomial UL**", "ApproxZ LL***", "ApproxZ UL***")
		} else { colnames(estResult) <- c("Variable", "Num Success", "Num Trial", "Proportion") }
	
		cat("Summary\n")
		printDataFrame(estResult)
		if (CI) {
			cat("* At ", confLevel*100,"% Confidence Level using Chi-Square Test\n", sep = "")
			cat("** At ", confLevel*100,"% Confidence Level using Exact Binomial Test\n", sep = "")
			cat("*** At ", confLevel*100,"% Confidence Level using Approximate Z Test\n", sep = "")
		}
		cat("\n")
		
		cat("Test on One Proportion\n")
		printDataFrame(summaryResult)
		cat("* alternative hypothesis = '", alternative,"'\n\n", sep = "")
		return(list(summaryEstimate = estResult, summaryStatistic = summaryResult))
	} ## end stmt -- if (procedure == "one")

	if (procedure == "independent") {
		estResult <- NULL
		estResult1 <- NULL
		summaryResult <- NULL
		#if (length(successLevel) != length(varX)) { stop("varX and successLevel should be of the same length.")}
		#if (length(testVal) != length(varX)) {
		#	if (length(testVal) < length(varX)) { testVal <- c(testVal, rep(0.5,(length(varX) - length(testVal))))
		#	} else { testVal <- testVal[1:length(varX)] }
		#}

		for (i in (1:length(varX))) {
			tempTable <- table(data[, varX[i]], data[,grp])
			if (nrow(tempTable) != 2) { 
                    cat(paste("ERROR: The variable '",varX[i],"' should be a dichotomous variable.\n\n", sep = "")) 
                    next
			} else { 
				if (any(is.na(match(rownames(tempTable), c("0","1"))))) {
					cat(paste("ERROR: The variable '",varX[i],"' should be a dichotomous variable with value 1 for success and 0 for failure.\n\n", sep = ""))
                         next
				} else { rownames(tempTable) <- c("Failure", "Success") }
			}
			successCnt <- tempTable["Success",]
			totalCnt <- colSums(tempTable)
			estProp <- successCnt/totalCnt

			chiResult <- suppressWarnings(prop.test(successCnt, totalCnt, conf.level = confLevel, correct = FALSE))
			pbar <- (successCnt[1] + successCnt[2])/(totalCnt[1] + totalCnt[2])
			approxZ <- (estProp[1] - estProp[2])/sqrt(pbar*(1-pbar)*((1/totalCnt[1])+(1/totalCnt[2])))

			if (alternative == "less") {	approxZ.pval <- pnorm(approxZ)
			} else {
				if (alternative == "greater") { approxZ.pval <- pnorm(approxZ, lower.tail = FALSE)
				} else { approxZ.pval <- 2*pnorm(-abs(approxZ)) }
			}

			if (i == 1) {
				if (CI) { 
					approxZWidth <- qnorm(1-((1-confLevel)/2))*sqrt(pbar*(1-pbar)*((1/totalCnt[1])+(1/totalCnt[2])))
					if (length(successCnt) == 2) {
						estResult <- data.frame(Variable = varX[i], paste("Diff(", paste(names(successCnt), collapse = "-"),")", sep = ""), 
									DiffEstimate = (estProp[1] - estProp[2]), chiResult$conf.int[1], UL1 = chiResult$conf.int[2],
									LL2 = (estProp[1] - estProp[2])- approxZWidth, UL2 = (estProp[1] - estProp[2]) + approxZWidth)
					} 
				} else { if (length(successCnt) == 2) { estResult <- data.frame(Variable = varX[i], paste("Diff(", paste(names(successCnt), collapse = "-"),")", sep = ""), DiffEstimate = (estProp[1] - estProp[2]))	}}
				estResult1 <- data.frame(Variable = varX[i], names(successCnt), t(rbind(successCnt, totalCnt, estProp)))
				summaryResult <- data.frame(Variable = varX[i],Procedure = "Chi-Square", Statistics = chiResult$statistic, Prob = chiResult$p.value)
			} else {
				if (CI) { 
					approxZWidth <- qnorm(1-((1-confLevel)/2))*sqrt(pbar*(1-pbar)*((1/totalCnt[1])+(1/totalCnt[2])))
					if (length(successCnt) == 2) {
						estResult <- rbind(estResult,
							 data.frame(Variable = varX[i], paste("Diff(", paste(names(successCnt), collapse = "-"),")", sep = ""), 
									DiffEstimate = (estProp[1] - estProp[2]), chiResult$conf.int[1], UL1 = chiResult$conf.int[2],
									LL2 = (estProp[1] - estProp[2])- approxZWidth, UL2 = (estProp[1] - estProp[2]) + approxZWidth))	
					}
					
				} else { if (length(successCnt) == 2) { estResult <- rbind(estResult, data.frame(Variable = varX[i], paste("Diff(", paste(names(successCnt), collapse = "-"),")", sep = ""), DiffEstimate = (estProp[1] - estProp[2])))	}}
				estResult1 <- rbind(estResult1, data.frame(Variable = varX[i], names(successCnt), t(rbind(successCnt, totalCnt, estProp))))
				summaryResult <- rbind(summaryResult,data.frame(Variable = varX[i],Procedure = "Chi-Square", Statistics = chiResult$statistic, Prob = chiResult$p.value))
			} 
			summaryResult <- rbind(summaryResult, data.frame(Variable = varX[i], Procedure = "Approx Z", Statistics = approxZ, Prob = approxZ.pval))
		}
		rownames(summaryResult) <- 1:nrow(summaryResult)
		rownames(estResult1) <- 1:nrow(estResult1)
		colnames(estResult1) <- c("Variable", grp, "Num Success", "Num Trial", "Proportion")
		if (length(successCnt) == 2) {
			if (CI) { colnames(estResult) <- c("Variable", grp, "Proportion Diff", "ChiSq LL*", "ChiSq UL*", "Approx Z LL**", "Approx Z UL**")
			} else { colnames(estResult) <- c("Variable", grp, "Proportion Diff") }
		}
		cat("Summary\n")
		printDataFrame(estResult1)
		cat("\n\n")
		if (length(successCnt) == 2) {
			printDataFrame(estResult)
			if (CI) {
				cat("* At ", confLevel*100,"% Confidence Level using Chi-Square Test\n", sep = "")
				cat("** At ", confLevel*100,"% Confidence Level using Approximate Z Test\n\n", sep = "")
			} else {
                    cat("\n")
			}
		}
		
		if (length(successCnt) == 2) { cat("Test on Two Independent Proportions\n")
		} else { cat("Test on Several Independent Proportions\n") }
		printDataFrame(summaryResult)
		#cat("* alternative = ", alternative,"\n\n")
		if (is.null(estResult)) { return(list(summaryEstimate = estResult1, summaryStatistic = summaryResult))
		} else { return(list(summaryEstimate = estResult1, summaryEstimateDiff = estResult, summaryStatistic = summaryResult)) }
		
	} ## end stmt -- if (procedure == "independent")

	if (procedure == "paired") {
		if (length(varX) != length(varY)) { stop("The arguments 'varX' and 'varY' should be of the same length") }
		summaryResult <- NULL
		for (i in (1:length(varX))) {
			if (length(unique(data[,varX[i]])) != 2 || length(unique(data[,varY[i]])) != 2) {
                    #cat("ERROR: The arguments 'varX' and 'varY' should both be dichotomous variables.")
                    cat(paste("ERROR: The variables '",varX[i],"' and '",varY[i],"' should both be dichotomous varibales.\n\n", sep = ""))
                    next
			}
			tempTable <- table(data[,c(varX[i], varY[i])])
			if (all(colnames(tempTable) == c("0","1"))) {
				colnames(tempTable) <- c("Failure", "Success")
			}
			if (all(rownames(tempTable) == c("0","1"))) {
				rownames(tempTable) <- c("Failure", "Success")
			}
			result <- mcnemar.test(tempTable, correct = FALSE)

			#cat(varX[i], "X", varY[i],"Contingency Table\n")
			cat("Contingency Table of", varX[i], "by", varY[i],"\n")
			printTable(tempTable)
			cat("\n\nTest on Two Related Proportion \nMcNemar'S Chi-Squared Test\n", sep = "")
			theWidth <- nchar("STATISTICS") + nchar(round(result$statistic[[1]],0))+ 10
			cat(paste(rep("-", theWidth), collapse = ""), "\n",sep = "")		
			cat(formatC("Statistics", width = 12, format = "s", flag = "-"), formatC(result$statistic[[1]], width = nchar(round(result$statistic[[1]],0))+ 8, digits = 4, format = "f"),"\n",sep = "")
			cat(formatC("DF", width = 12, format = "s", flag = "-"), formatC(result$parameter[[1]], width = nchar(round(result$statistic[[1]],0))+ 8, format = "d"),"\n", sep = "")
			cat(formatC("Prob", width = 12, format = "s", flag = "-"), formatC(result$p.value[[1]], width = nchar(round(result$statistic[[1]],0))+ 8, digits = 4, format = "f"),"\n",sep = "")
			cat(paste(rep("-", theWidth), collapse = ""),sep = "")		
			cat("\n\n")
		}
	}
}