# ---------------------------------------------------------------------------------------------
# R-CropStat Beta Version: Function for ANALYZE - REGRESSION AND CORRELATION SUBMENU
# ---------------------------------------------------------------------------------------------
# BivariateCorrelationTest: Perform correlation analysis using pearson, spearman, and
# Created by: Alaine A. Gulles 03.01.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.19.2011
# ---------------------------------------------------------------------------------------------
# ARGUMENTS: x = numeric data matrix
#		 method = a character string indicating which correlation 
#                     coefficient is to be used for the test. User can specify
#			    all three correlation analysis. By default, performs the
#			    pearson correlation analysis.
#		 alternative = indicates the alternative hypothesis to be used and
#                          must be one of the following "two.sided", "greater"
#                          or "less".
# ---------------------------------------------------------------------------------------------

BivariateCorrelationTest <- function(data, var, method = "pearson", alternative = "two.sided", statistics = FALSE) UseMethod("BivariateCorrelationTest")

BivariateCorrelationTest.default <- function(data, var, method = "pearson", alternative = "two.sided", statistics = FALSE) {
	if (is.character(data)) {
		nameData <- data
		if (!exists(data)) { stop(paste("The object '", nameData,"' not found.", sep = "")) }
		tempData <- eval(parse(text = data))
	} else {
		nameData <- paste(deparse(substitute(data)))
		tempData <- data
	}
	if (!is.data.frame(tempData)) { stop("The object 'data' should be of type data frame.") }
	if (!is.character(var)) { stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }

	procedure = c("pearson", "spearman", "kendall")
	method <- procedure[match(method, procedure)]
	if (all(is.na(method))) stop("method does not match")
	alternative <- match.arg(alternative, c("two.sided",  "greater", "less"))

	pairs.var <- t(combn(var,2))
	result <- list()
	options(width = 5000, digits = 6, warn = -1)
	if (statistics) { DescriptiveStatistics(data, var, grp = NULL, statistics = c("nnmiss", "min", "max", "mean", "sd")); cat("\n\n")	}

	if (alternative == "two.sided") { sub.title <- "Prob > |r|"
	} else { 
		if (alternative == "greater") { sub.title <- "Prob > r" 
		} else { sub.title <- "Prob < r" }
	}
	cat("CORRELATION ANALYSIS","\n\n")

	for (i in (1:length(method))){
		pval.matrix <- as.table(matrix(NA, nrow = length(var), ncol = length(var), dimnames = list(var, var)))
		est.matrix <- as.table(matrix(1, nrow = length(var), ncol = length(var), dimnames = list(var, var)))
		nobs.matrix <- est.matrix
		for (j in (1:length(var))) { nobs.matrix[j,j] <- length(na.omit(tempData[,var[j]])) }
		for (j in (1:nrow(pairs.var))){
			theData <- na.omit(tempData[,c(pairs.var[j,1], pairs.var[j,2])])
			result.corr <- cor.test(theData[,1], theData[,2], method = method[i], alternative = alternative, exact = FALSE)
			est.matrix[match(pairs.var[j,1], row.names(est.matrix)),match(pairs.var[j,2], row.names(est.matrix))] <- result.corr$estimate
			est.matrix[match(pairs.var[j,2], row.names(est.matrix)),match(pairs.var[j,1], row.names(est.matrix))] <- result.corr$estimate
			pval.matrix[match(pairs.var[j,1], row.names(pval.matrix)),match(pairs.var[j,2], row.names(pval.matrix))] <- result.corr$p.value
			pval.matrix[match(pairs.var[j,2], row.names(pval.matrix)),match(pairs.var[j,1], row.names(pval.matrix))] <- result.corr$p.value
			nobs.matrix[match(pairs.var[j,1], row.names(nobs.matrix)),match(pairs.var[j,2], row.names(nobs.matrix))] <- nrow(theData)
			nobs.matrix[match(pairs.var[j,2], row.names(nobs.matrix)),match(pairs.var[j,1], row.names(nobs.matrix))] <- nrow(theData)
		}
		cat(result.corr$method, ", ", sub.title,"\n\n", sep = "") 
		print(printCorrMatrix(est.matrix, pval.matrix, nobs.matrix), row.names = TRUE, right = TRUE)
		cat("\n\n")
		result[[i]] <- list(procedure = result.corr$method, estimate = est.matrix, pvalue = pval.matrix, nobs = nobs.matrix, alternative = alternative)
	}
	return(invisible(result))
} ### end stmt -- bivariateCorrelationTest
