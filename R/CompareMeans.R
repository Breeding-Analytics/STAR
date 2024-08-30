# -------------------------------------------------------------------------------------
# R-CropStat Beta Version: 
# -------------------------------------------------------------------------------------
# CompareMeans
# Created by: Alaine A. Gulles 09.30.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.04.2012
# -------------------------------------------------------------------------------------

CompareMeans <- function(
	data, varX, varY = NULL, grp = NULL, testVal = 0,
	procedure = c("one","paired","independent"), 
	alternative = c("two.sided", "greater", "less"), statistics = FALSE,
	CI = FALSE, confLevel = 0.95, normality = NULL, alpha = 0.05) UseMethod("CompareMeans")

CompareMeans.default <- function(
	data, varX, varY = NULL, grp = NULL, testVal = 0,
	procedure = c("one","paired","independent"), 
	alternative = c("two.sided", "greater", "less"), statistics = FALSE,
	CI = FALSE, confLevel = 0.95, normality = NULL, alpha = 0.05) {

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
	if (!is.data.frame(tempData)) { stop("The argument 'data' should be of type data frame.") }
	if (!is.character(varX)) 	{ stop(paste("The argument 'varX' should be a character vector.", sep = "")) }
	if (any(is.na(match(varX, names(tempData))))) { stop("At least one item in the character vector 'varX' does not match any variable name in the dataset.") }

	if (procedure == "one") {
		if (!missing(testVal) && (length(testVal) != 1 || is.na(testVal))) { stop("The argument 'testVal' must be a single number.") }
	}
   	if (!missing(confLevel) && (length(confLevel) != 1 || !is.finite(confLevel) || confLevel < 0 || confLevel > 1))  stop("The argument 'confLevel' must be a single number between 0 and 1")
	alternative <- match.arg(alternative)
	procedure <- match.arg(procedure)
	if (!is.null(normality)) {
		methodNormality = c("swilk", "sfrancia", "ks", "cramer", "anderson")
		normality <- methodNormality[match(normality, methodNormality)]
		if (all(is.na(normality))) normality <- NULL
	}

	ttestResult <- NULL
	statisticsResult <- NULL
	CIResult <- NULL
	if (alternative == "two.sided") { pvalueLabel <- "Pr(>|t|)"	} else {
		if (alternative == "greater") pvalueLabel <- "Pr(> t)" else pvalueLabel <- "Pr(< t)"
	}
	
	if (procedure == "one") {
		for (i in (1:length(varX))) {
			result <- t.test(tempData[,varX[i]], mu = testVal, alternative = alternative, conf.level = confLevel)
			ttestResult <- rbind(ttestResult, data.frame(Variable = varX[i], DF = result[[2]], TESTSTAT = result[[1]], PVALUE = result[[3]]))
			if (CI) { 
				tempResult <- t.test(tempData[,varX[i]], mu = testVal, alternative = "two.sided", conf.level = confLevel)
				CIResult <- rbind(CIResult, data.frame(Variable = varX[i], Lower_CI = tempResult[[4]][1], Mean = result$estimate, Upper_CI = tempResult[[4]][2])) 
			}
		}
		rownames(ttestResult) <- 1:nrow(ttestResult)
		colnames(ttestResult) <- c("Variable", "DF", "t Value", pvalueLabel)

		capture.output(statisticsResult <- DescriptiveStatistics(data = tempData, var = varX, grp = NULL, statistics = c("nnmiss", "mean", "sd", "se.mean")))
		if (CI) {
			statisticsResult <- data.frame(statisticsResult[,1:2], CIResult[,2:4], statisticsResult[,4:5])
			colnames(CIResult) <- c("Variable", "Lower CI*", "Mean", "Upper CI*")
			colnames(statisticsResult)[c(2,3,5)] <- c("N", "Lower CI*", "Upper CI*")
		}
		if (!is.null(normality)) {
			NormalityTest(data = tempData, var = varX, grp = NULL, method = normality) 
			cat("\n")
		}
		if (statistics) {
			cat("Descriptive Statistics\n")
			printDataFrame(statisticsResult)
			if (CI) { cat("* At ", confLevel*100,"% Confidence Level.\n\n\n", sep = "") } else { cat("\n\n") }
		} else {
			if (CI) {
				cat("Confidence Interval of the Mean\n")
				printDataFrame(CIResult)
				cat("* At ", confLevel*100,"% Confidence Level.\n\n", sep = "")
			}
		}

		cat("One Sample t-Test, h0: mean = ", testVal,"\n", sep = "")
		printDataFrame(ttestResult)
		
		if (CI) {
			return(invisible(list(statistics = ttestResult, 
					estimate = statisticsResult, 
					confLevel = confLevel,
					nullValue = testVal,
					alternative = alternative, 
					method = "One Sample t-Test")))
		} else {
			return(invisible(list(statistics = ttestResult, 
					estimate = statisticsResult, 
					nullValue = testVal,
					alternative = alternative, 
					method = "One Sample t-Test")))
		}
	}


	if (procedure == "paired") {
		if (is.null(varY)) { stop("The argument 'varY' is required.") }
		if (!is.character(varY)) 	{ stop(paste("The argument 'varY' should be a character vector.", sep = "")) }
		if (any(is.na(match(varY, names(tempData))))) { stop("At least one item in the character vector 'varY' does not match any variable name in the dataset.") }
		if (length(varX) != length(varY)) { stop("The arguments 'varX' and 'varY' should be of the same length.") }
		varDiff = data.frame(tempData[,varX] - tempData[,varY])
		names(varDiff) <- varX
		varLabel <- trimStrings(strsplit(paste(varX, collapse = " * ", varY, sep = " - "), split = "\\*")[[1]])

		for (i in (1:ncol(varDiff))) {
			result <- t.test(x = tempData[,varX[i]], y = tempData[,varY[i]], mu = testVal, alternative = alternative, conf.level = confLevel, paired = TRUE)
			ttestResult <- rbind(ttestResult, data.frame(Difference = varLabel[i], DF = result[[2]], TESTSTAT = result[[1]], PVALUE = result[[3]]))
			if (CI) { 
				tempResult <- t.test(x = tempData[,varX[i]], y = tempData[,varY[i]], mu = testVal, alternative = "two.sided", conf.level = confLevel, paired = TRUE)
				CIResult <- rbind(CIResult, data.frame(Difference = varLabel[i], Lower_CI = tempResult[[4]][1], Mean = result[[5]], Upper_CI = tempResult[[4]][2])) 
			}
		}
		rownames(ttestResult) <- 1:nrow(ttestResult)
		colnames(ttestResult)[3:4] <- c("t Value", pvalueLabel)

		capture.output(statisticsResult <- DescriptiveStatistics(data = varDiff, var = varX, grp = NULL, statistics = c("nnmiss", "mean", "sd", "se.mean")))
		colnames(statisticsResult)[1] <- "Difference"
		statisticsResult[,"Difference"] <- varLabel

		if (CI) {
			statisticsResult <- data.frame(statisticsResult[1:2], CIResult[2:4], statisticsResult[4:5])
			colnames(CIResult) <- c("Difference", "Lower CI*", "Mean", "Upper CI*")
			colnames(statisticsResult)[c(2,3,5)] <- c("N", "Lower CI*", "Upper CI*")
		}

		if (!is.null(normality)) {
			capture.output(resultNormality <- NormalityTest(data = varDiff, var = varX, grp = NULL, method = normality))
			colnames(resultNormality)[1] <- "Difference"
			resultNormality[,"Difference"] <- varLabel
			cat("Test for Normality\n")
			printDataFrame(resultNormality)
			cat("\n")
		}

		if (statistics) {
			cat("Descriptive Statistics\n")
			printDataFrame(statisticsResult)
			if (CI) {  cat("* At ", confLevel*100,"% Confidence Level.\n\n", sep = "") }
		} else {
			if (CI) {
				cat("Confidence Interval of the Mean Difference\n")
				printDataFrame(CIResult)
				cat("* At ", confLevel*100,"% Confidence Level.\n\n", sep = "")
			}
		}
		cat("\n")
		cat("Paired Sample t-Test, h0: mean diff = ", testVal,"\n", sep = "")
		printDataFrame(ttestResult)
		cat("\n")

		if (CI) { 
			return(invisible(list(statistics = ttestResult, 
					estimate = statisticsResult, 
					confLevel = confLevel,
					nullValue = testVal,
					alternative = alternative, 
					method = "Paired Sample t-Test"))) 
		} else {
			return(invisible(list(statistics = ttestResult, 
					estimate = statisticsResult, 
					nullValue = testVal,
					alternative = alternative, 
					method = "Paired Sample t-Test"))) 
		}
	}

	if (procedure == "independent") {
		if (is.null(grp)) { stop("The argument 'grp' is required.") }
		if (!is.character(grp)) 	{ stop(paste("The argument 'grp' should be a character string.", sep = "")) }
		if (length(grp) != 1) { stop("The argument 'grp' should be a single character string.") }
		if (any(is.na(match(grp, names(tempData))))) { stop("The argument 'grp' does not match any variable name in the dataset.") }
          tempData[,grp] <- factor(tempData[,grp])
		if (nlevels(tempData[,grp]) != 2) { stop("The argument 'grp' has more than two levels.") }
	   	if (!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) || alpha < 0 || alpha > 1))  stop("The argument 'alpha' must be a single number between 0 and 1")

		equalVarResult <- NULL
		for (i in (1:length(varX))) {
			capture.output(statisticsTemp <- DescriptiveStatistics(data = tempData, var = varX[i], grp = grp, statistics = c("nnmiss", "mean", "sd", "se.mean")))
			statisticsTemp[,2] <- as.character(statisticsTemp[,2])
			if (CI) {
				ttab <- qt(1-((1-confLevel)/2), statisticsTemp[,"N_NonMissObs"]-1)
				LL <- statisticsTemp[,"Mean"] - (ttab*(statisticsTemp[1,"StdDev"]/sqrt(statisticsTemp[,"N_NonMissObs"])))
				UL <- statisticsTemp[,"Mean"] + (ttab*(statisticsTemp[1,"StdDev"]/sqrt(statisticsTemp[,"N_NonMissObs"])))
				statisticsTemp <- data.frame(statisticsTemp[,1:3], LowerCI = LL, statisticsTemp[4], UpperCI = UL, statisticsTemp[,5:6])
			}
			statisticsTemp <- rbind(statisticsTemp, statisticsTemp[1,])
			foldedResult <- FoldedFTest(tempData[,varX[i]], tempData[,grp])
			if (foldedResult$pvalue < alpha) {
				equalVar <- FALSE 
				ttestMethod <- "Satterthwaite"
				varMethod <- "Unequal"
			} else { 
				equalVar <- TRUE
				ttestMethod <- "Pooled"
				varMethod <- "Equal"
			}
			equalVarResult <- rbind(equalVarResult, data.frame(Variable = varX[i], Method = foldedResult[[4]], NumDF = foldedResult[[2]][[1]], DenDF = foldedResult[[2]][[2]], FValue = foldedResult[[1]], PValue = foldedResult[[3]]))
			result <- t.test(formula(paste(varX[i],"~", grp)), data = tempData, mu = testVal, alternative = alternative, conf.level = confLevel, var.equal = equalVar)
			ttestResult <- rbind(ttestResult, data.frame(Variable = varX[i], Method = ttestMethod, Variances = varMethod,  DF = result$parameter[[1]], tValue = result$statistic, PValue = result$p.value))
			statisticsTemp[3,2] <- paste("Diff(", paste(levels(tempData[,grp]), collapse = "-"),")", sep = "")
			statisticsTemp[3,3] <- NA
			statisticsTemp[3,"Mean"] <- result$estimate[[1]] - result$estimate[[2]]
			pooledVar <- (((statisticsTemp[1,"N_NonMissObs"]-1) * statisticsTemp[1,"StdDev"]**2) + ((statisticsTemp[2,"N_NonMissObs"]-1) * statisticsTemp[2,"StdDev"]**2))/(sum(statisticsTemp[,"N_NonMissObs"], na.rm = TRUE) - 2)
			statisticsTemp[3,"StdDev"] <- sqrt(pooledVar)
			statisticsTemp[3,"SE_Mean"] <- sqrt(pooledVar*((1/statisticsTemp[1,"N_NonMissObs"])+(1/statisticsTemp[2,"N_NonMissObs"])))

			if (CI) {
				tempResult <- t.test(formula(paste(varX[i],"~", grp)), data = tempData, mu = testVal, alternative = "two.sided", conf.level = confLevel, var.equal = equalVar)
				statisticsTemp[3, "LowerCI"] <- tempResult$conf.int[1]
				statisticsTemp[3, "UpperCI"] <- tempResult$conf.int[2]
			}
			statisticsResult <- rbind(statisticsResult, statisticsTemp)
		}
		if (CI) { colnames(statisticsResult)[c(3,4,6)] <- c("N", "Lower CI*", "Upper CI*")
		} else { colnames(statisticsResult)[c(3)] <- c("N") }
		
		colnames(equalVarResult) <- c("Variable", "Method", "Num DF", "Den DF", "F Value", "Pr(> F)")
		#equalVarResult[,c(3:4)] <- as.integer(equalVarResult[,c(3:4)])
		colnames(ttestResult) <- c("Variable", "Method*", "Variances", "DF", "t Value", pvalueLabel)
		rownames(ttestResult) <- 1:nrow(ttestResult)

		if (!is.null(normality)) {
			NormalityTest(data = tempData, var = varX, grp = grp, method = normality)
			cat("\n")
		}

		if (statistics) {
			cat("Descriptive Statistics\n")
			printDataFrame(statisticsResult)
			if (CI) { cat("* At ", confLevel*100,"% Confidence Level.\n\n", sep = "") } else { cat("\n\n") }
		} else {
			if (CI) {
				cat(toupper("Confidence Interval of the Mean\n"))
				temp <- statisticsResult[,c(1,2,4:6)]
				printDataFrame(temp)
				cat("* At ", confLevel*100,"% Confidence Level.\n\n", sep = "")
			} else { cat("\n\n") }
		}
		cat("Homogeneity of Variances\n")
		printDataFrame(equalVarResult)
		cat("\n")		
		
		cat("Two Independent Sample t-Test, h0: mean diff = ", testVal,"\n", sep = "")
		printDataFrame(ttestResult)
		cat("* At ", alpha," level of significance.\n\n", sep = "")		

		if (CI) {
			return(invisible(list(statistics = ttestResult, 
					equalVariance = equalVarResult,
					alpha = alpha,
					estimate = statisticsResult, 
					confLevel = confLevel,
					nullValue = testVal,
					alternative = alternative, 
					method = "Two Independent Sample t-Test")))

		} else {
			return(invisible(list(statistics = ttestResult, 
					equalVariance = equalVarResult,
					alpha = alpha,
					estimate = statisticsResult, 
					nullValue = testVal,
					alternative = alternative, 
					method = "Two Independent Sample t-Test")))

		}
	}
} 