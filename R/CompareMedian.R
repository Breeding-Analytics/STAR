# -------------------------------------------------------------------------------------
# R-CropStat Beta Version: 
# -------------------------------------------------------------------------------------
# CompareMedians
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.04.2012
# -------------------------------------------------------------------------------------

CompareMedian <- function(data, varX, varY = NULL, grp = NULL, procedure = c("one", "paired", "independent"), altHypo = c("two.sided", "less", "greater"), hypoValue = 0, confInt = FALSE, confLevel = 0.95) UseMethod("CompareMedian")

CompareMedian.default <- function(data, varX, varY = NULL, grp = NULL, procedure = c("one", "paired", "independent"), altHypo = c("two.sided", "less", "greater"), hypoValue = 0, confInt = FALSE, confLevel = 0.95) {
	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else { nameData <- paste(deparse(substitute(data))) }

	altHypo <- match.arg(altHypo)
	procedure <- match.arg(procedure)
	tempTable <- NULL

	# perform wilcoxon signed rank test
	if (procedure == "one") {
		for (i in (1:length(varX))) {
			tempResult <- wilcox.test(x = data[,varX[i]], alternative = altHypo, mu = hypoValue, paired = FALSE, conf.int = TRUE, conf.level = confLevel)
			if (confInt) {
				tempTable <- rbind(tempTable, 
							 data.frame(Variable = varX[i], Statistics = tempResult$statistic,
									Prob = tempResult$p.value, LL_CI = tempResult$conf.int[1], 
			     						Median = tempResult$estimate, UL_CI = tempResult$conf.int[2]))
			} else {
				tempTable <- rbind(tempTable, 
						 	 data.frame(Variable = varX[i], Statistics = tempResult$statistic[[1]],
							 Prob = tempResult$p.value, Median = tempResult$estimate))
			}
		}
		if (altHypo == "two.sided") { colnames(tempTable)[3] <- paste("Prob > |",attr(tempResult$statistic, "names"),"|", sep = "") 
		} else {	if (altHypo == "greater") { colnames(tempTable)[3] <- paste("Prob > ",attr(tempResult$statistic, "names"), sep = "") 
				} else { colnames(tempTable)[3] <- paste("Prob < ",attr(tempResult$statistic, "names"), sep = "") }
		}
	
		if (confInt) { 
			colnames(tempTable)[4] <- paste(colnames(tempTable)[4], "*") 
			colnames(tempTable)[6] <- paste(colnames(tempTable)[6], "*") 
		}
	
		cat("Sigh Test: Md = ", hypoValue, "\n", sep = "")
		printDataFrame(tempTable)
		if (confInt) {
			cat("* ", confLevel*100, "% Confidence Interval\n\n", sep = "")
		}
		return(invisible(tempTable))
	} ## end stmt -- if (procedure == "one")

	if (procedure == "paired") {
	     if (is.null(varY)) { stop("The argument 'varY' is required.") }
	     if (length(varX) != length(varY)) { stop("The arguments 'varX' and 'varY' should be of the same length.") }
	     
		for (i in (1:length(varX))) {
			tempResult <- wilcox.test(x = data[,varX[i]], y = data[,varY[i]], alternative = altHypo, mu = hypoValue, paired = TRUE, conf.int = TRUE, conf.level = confLevel)
			if (confInt) {
				tempTable <- rbind(tempTable, 
						 data.frame(Variable = paste(varX[i], "-", varY[i]), Statistics = tempResult$statistic,
								Prob = tempResult$p.value, LL_CI = tempResult$conf.int[1], 
		     					Median = tempResult$estimate, UL_CI = tempResult$conf.int[2]))
			} else {
				tempTable <- rbind(tempTable, 
						 data.frame(Variable = paste(varX[i], "-", varY[i]), Statistics = tempResult$statistic,
								Prob = tempResult$p.value, Median = tempResult$estimate))
			}
		} 

		if (altHypo == "two.sided") { colnames(tempTable)[3] <- paste("Prob > |",attr(tempResult$statistic, "names"),"|", sep = "") 
		} else {	if (altHypo == "greater") { colnames(tempTable)[3] <- paste("Prob > ",attr(tempResult$statistic, "names"), sep = "") 
				} else { colnames(tempTable)[3] <- paste("Prob < ",attr(tempResult$statistic, "names"), sep = "") }
		}

		if (confInt) { 
			colnames(tempTable)[4] <- paste(colnames(tempTable)[4], "*") 
			colnames(tempTable)[6] <- paste(colnames(tempTable)[6], "*") 
		}
	
		cat("Wilcoxon Signed Rank Test\n\n")
		printDataFrame(tempTable)
		if (confInt) {
			cat("* ", confLevel*100, "% Confidence Interval\n\n", sep = "")
		}
		return(invisible(tempTable))
	} ## end stmt -- if (procedure == "paired")

	if (procedure == "independent") {
		newGrp <- NULL
		for (i in (1:length(grp))) {
			data[,grp[i]] <- factor(data[,grp[i]])
			if (nlevels(data[,grp[i]]) == 2) {
				newGrp <- c(newGrp, grp[i])
			}
			
		}
		if (is.null(newGrp)) {
			stop("The argument 'grp' should only have two levels for Mann-Whitney Test to be performed.")
		} else {
			grp <- newGrp
		}
		for (i in (1:length(grp))) {
			for (j in (1:length(varX))) {
				tempResult <- wilcox.test(data[,varX[j]]~ data[,grp[i]], alternative = altHypo, mu = hypoValue, conf.int = TRUE, conf.level = confLevel)
				if (confInt) {
					tempTable <- rbind(tempTable, 
							 data.frame(Variable = varX[j], Statistics = tempResult$statistic,
									Prob = tempResult$p.value, LL_CI = tempResult$conf.int[1], 
				     					Median = tempResult$estimate, UL_CI = tempResult$conf.int[2]))
				} else {
					tempTable <- rbind(tempTable, 
								 data.frame(Variable = varX[j], Statistics = tempResult$statistic,
								 Prob = tempResult$p.value, Median = tempResult$estimate))
				}
			}
		} 


		if (altHypo == "two.sided") { colnames(tempTable)[3] <- paste("Prob > |",attr(tempResult$statistic, "names"),"|", sep = "") 
		} else {	if (altHypo == "greater") { colnames(tempTable)[3] <- paste("Prob > ",attr(tempResult$statistic, "names"), sep = "") 
			} else { colnames(tempTable)[3] <- paste("Prob < ",attr(tempResult$statistic, "names"), sep = "") }
		}

		if (confInt) { 
			colnames(tempTable)[4] <- paste(colnames(tempTable)[4], "*") 
			colnames(tempTable)[6] <- paste(colnames(tempTable)[6], "*") 
		}
	
		cat("Mann-Whitney (Wilcoxon Rank Sum) Test\n\n")
		printDataFrame(tempTable)
		if (confInt) {
			cat("* ", confLevel*100, "% Confidence Interval\n\n", sep = "")
		}
		return(invisible(tempTable))
	} ## end stmt -- if (procedure == "independent")

}