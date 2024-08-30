# -----------------------------------------------------------------------
# LinearRegressionAnalysis: GUI and function for Linear Regression Analysis
# Created by: Alaine A. Gulles 10.21.2010
# Modified by: Alaine A. Gulles 09.24.2012 
# -----------------------------------------------------------------------

LinearRegressionAnalysis <- function(data, depVar, indepVar, constant = TRUE, 
                                     statistics = TRUE, confInt = FALSE, confLevel = 0.95, 
                                     covMatrix = FALSE, normality = NULL, 
                                     heteroskedasticity = NULL, autoCorr = FALSE, 
                                     VIF = FALSE, COOKS = FALSE, leverage = FALSE) UseMethod("LinearRegressionAnalysis")
     
LinearRegressionAnalysis.default <- function(data, depVar, indepVar, constant = TRUE, 
                                             statistics = TRUE, confInt = FALSE, confLevel = 0.95, 
                                             covMatrix = FALSE, normality = NULL, 
                                             heteroskedasticity = NULL, autoCorr = FALSE, 
                                             VIF = FALSE, COOKS = FALSE, leverage = FALSE) {

	if (is.character(data)) {
		nameData <- data
		data <- eval(parse(text = data))
	} else { nameData <- paste(deparse(substitute(data))) }

	if (!is.null(normality)) { 
		availMethod <- c("swilk", "sfrancia", "ks", "cramer", "anderson")
		normality <- availMethod[na.omit(match(normality, availMethod))]
		if (length(normality) == 0) { normality <- availMethod[1] }
	}

	if (!is.null(heteroskedasticity)) { 
		availMethod <- c("bpagan", "gquandt")
		heteroskedasticity <- availMethod[na.omit(match(heteroskedasticity, availMethod))]
		if (length(heteroskedasticity) == 0) { heteroskedasticity <- NULL }
	}

	if (!is.numeric(confLevel)) { confLevel <- 0.95 } else { if (confLevel >= 1 || confLevel < 0.9) confLevel <- 0.95 }

	tempData <- list()
	myfit <- list()
	diagPlot <- list()
	numVar <- NULL
	factorVar <- NULL
     rownames(data) <- 1:nrow(data)
     	
     if (statistics) {
      	for (i in (1:length(indepVar))) {
            	if (is.numeric(data[,indepVar[i]])) { numVar <- c(numVar, indepVar[i]) } else { factorVar <- c(factorVar, indepVar[i]) }     
     	}
          
     	DescriptiveStatistics(data = data, var = c(depVar, numVar), statistics = c("nnmiss", "min", "max", "mean", "sd", "se.mean")) 
     	cat("\n")
          
     	if (!is.null(factorVar)) {
          	cat("Observed Frequency for Qualitative Variables\n")
          	for (i in (1:length(factorVar))) {
               	tempTable <- data.frame(table(data[,factorVar[i]]))
               	newTable <- data.frame("Freq", t(tempTable[,ncol(tempTable)]))
               	colnames(newTable) <- c(factorVar[i], as.character(t(tempTable[,1])))
               	printDataFrame(newTable)
               	cat("\n")
          	}
     	}
     }

	for (i in (1:length(depVar))) {
		tempResult <- NULL
		tempAOV <- NULL
		modelSummary <- NULL
		residData <- NULL

         	if (constant) { theModel <- paste(depVar[i], "~", paste(indepVar, collapse = " + "))
         	} else { theModel <- paste(depVar[i], "~", paste(indepVar, collapse = " + "), "- 1") }
		
		command <- paste("myfit[[",i,"]] <- lm(",theModel,", data = ",nameData,", y = TRUE)", sep = "")
		eval(parse(text = command))
		tempResult <- summary(myfit[[i]])
		tempAOV <- anova(myfit[[i]])

		cat(toupper("Linear Regression Analysis"), "\n")
		cat("Model Fitted: ", theModel,"\n\n")
		#cat("Dependent Variable: ", depVar[i],"\n\n")
		
		# DISPLAY THE MODEL FIT
		cat("Analysis of Variance Table\n")	
		#temp <- rbind(tempAOV[1,],tempAOV)
		#temp[1,] <- cbind(sum(tempAOV[1:(nrow(tempAOV)-1),1]),
		#			sum(tempAOV[1:(nrow(tempAOV)-1),2]),
		#			sum(tempAOV[1:(nrow(tempAOV)-1),2])/tempResult$fstatistic[[2]],
		#			tempResult$fstatistic[[1]],
		#			pf(tempResult$fstatistic[[1]],tempResult$fstatistic[[2]],tempResult$fstatistic[[3]], lower.tail = FALSE))
          #tempRNames <- paste(" ",rownames(tempAOV)[1:(nrow(tempAOV)-1)])
          #rownames(temp) <- c("Model", tempRNames, rownames(tempAOV)[nrow(tempAOV)])

		temp <- tempAOV[c(1, nrow(tempAOV)),]
		temp[1,] <- cbind(sum(tempAOV[1:(nrow(tempAOV)-1),1]), 
                            sum(tempAOV[1:(nrow(tempAOV)-1),2]), 
                            sum(tempAOV[1:(nrow(tempAOV)-1),2])/tempResult$fstatistic[[2]],
                            tempResult$fstatistic[[1]],
                            pf(tempResult$fstatistic[[1]],tempResult$fstatistic[[2]],tempResult$fstatistic[[3]], lower.tail = FALSE))
		rownames(temp) <- c("Model", rownames(tempAOV)[nrow(tempAOV)])
		printAOVTable(ConstructAOVTable(temp))
		cat("\n")
			
		modelSummary <- data.frame(tempResult$sigma, mean(myfit[[i]]$y, na.rm = TRUE), cv(myfit[[i]]$y, na.rm = TRUE), tempResult$r.squared, tempResult$adj.r.squared)
		colnames(modelSummary) <- c("Root MSE", paste(depVar[i], "Mean"), "CV(%)", "R-Square", "Adj R-Sq")			
		cat("Model Summary:\n")
		printDataFrame(modelSummary)
		cat("\n")
			
		# DISPLAY THE PARAMETER ESTIMATE
		cat("Parameter Estimates:\n")
          if (constant) { rownames(tempResult$coefficient)[1] <- "Intercept" }
		if (confInt) {
		     #coefTable <- data.frame(rownames(tempResult$coefficient),tempResult$coefficient[,1:3], confint(myfit[[i]], level = confLevel), tempResult$coefficient[,4])     
               if (nrow(tempResult$coefficient) != 1) {
                    coefTable <- data.frame(rownames(tempResult$coefficient),tempResult$coefficient[,1:3], confint(myfit[[i]], level = confLevel), tempResult$coefficient[,4])     
               } else {
                    coefTable <- data.frame(rownames(tempResult$coefficient),t(tempResult$coefficient[,1:3]), confint(myfit[[i]], level = confLevel), tempResult$coefficient[,4])     
               }
			colnames(coefTable) <- c("Variable", attributes(tempResult$coefficient)$dimnames[[2]][1:3], "LL CI*","UL CI*", attributes(tempResult$coefficient)$dimnames[[2]][4])
               #names(coefTable)[ncol(coefTable)] <- "Pr > |t|"
			printDataFrame(coefTable)
			cat("* At ",(confLevel)*100,"% Confidence Interval\n\n", sep = "")
		} else {
			coefTable <- data.frame(rownames(tempResult$coefficient),tempResult$coefficient)
			colnames(coefTable) <- c("Variable", attributes(tempResult$coefficient)$dimnames[[2]])
			#names(coefTable)[ncol(coefTable)] <- "Pr > |t|"
			printDataFrame(coefTable)
			cat("\n")
		}

		if (covMatrix) {
			cat("Coefficient Variance-Covariance Matrix:\n")
			print(vcov(myfit[[i]])[2:nrow(vcov(myfit[[i]])),2:ncol(vcov(myfit[[i]]))], digits = 5)
			cat("\n\n")
		}

		if (!is.null(normality)) {
			residData <- data.frame(myfit[[i]]$residuals)
			colnames(residData) <- "residual"
			NormalityTest(data = residData, var = "residual", grp = NULL, method = normality)
			cat("\n\n")
		}

		if (!is.null(heteroskedasticity)) {
			heteroTable <- NULL
			if (!is.na(match("bpagan", heteroskedasticity))) {
				result <- bptest(myfit[[i]], studentize = FALSE)
				heteroTable <- data.frame(Method = result[[3]], DF = result[[2]], Statistic = names(result[[1]]), Value = result[[1]], pvalue = result[[4]])
			} 
			if (!is.na(match("gquandt", heteroskedasticity))) {
				result <- gqtest(myfit[[i]], alternative = "two.sided")
				if (!is.null(heteroTable)) {
					heteroTable <- data.frame(Method = heteroTable[,1], DF1 = heteroTable[,2],DF2 = NA, heteroTable[,3:ncol(heteroTable)])
					heteroTable <- rbind(heteroTable, data.frame(Method = result[[3]], DF1 = result[[2]][[1]], DF2 = result[[2]][[2]], Statistic = names(result[[1]]), Value = result[[1]], pvalue = result[[4]]))
				} else {
					heteroTable <- data.frame(Method = result[[3]], DF1 = result[[2]][[1]], DF2 = result[[2]][[2]], Statistic = names(result[[1]]), Value = result[[1]], pvalue = result[[4]])
				}
			}
			
               if (length(heteroskedasticity) == 1) {
                    pvalIndex <- match("pvalue", names(heteroTable))
                    valIndex <- match("Value", names(heteroTable))
                    names(heteroTable)[valIndex] <- paste(heteroTable[,"Statistic"], "Value")
                    names(heteroTable)[pvalIndex] <- paste("Pr(>",heteroTable[,"Statistic"], ")", sep = "")
                    heteroTable[,"Statistic"] <- NULL
                    cat("Test for Heteroskedasticity\n")
               } else {
                    pvalIndex <- match("pvalue", names(heteroTable))
                    pval <- heteroTable[,pvalIndex]
                    heteroTable <- data.frame(heteroTable[,1:(pvalIndex-1)], paste("Pr(>",heteroTable[,"Statistic"],")",sep = ""), heteroTable[,pvalIndex])
                    names(heteroTable)[pvalIndex:(pvalIndex+1)] <- c("Prob", "p Value")
                    cat("Tests for Heteroskedasticity\n")
               }
			printDataFrame(heteroTable)
			cat("\n\n")
		}
		
		if (autoCorr) {
			cat("Durbin Watson Test for Autocorrelation:\n")
			#print(durbinWatsonTest(myfit[[i]], alternative = "two.sided"))
			resultdw <- durbinWatsonTest(myfit[[i]], alternative = "two.sided")
			dwTable <- data.frame(Autocorrelation = resultdw[[1]], Statistic = resultdw[[2]], Prob = resultdw[[3]])
               names(dwTable)[3] <- "Pr(>|D-W|)"
               printDataFrame(dwTable)
			cat("\n\n")
		}

		if (VIF) {
			cat("Variance Inflation Factor:\n")
               	if (is.null(nrow(vif(myfit[[i]])))) { printDataFrame(data.frame(t(suppressWarnings(vif(myfit[[i]])))))
               	} else {
                    	vifTable <- data.frame(Variables = indepVar, suppressWarnings(vif(myfit[[i]])))
                    	colnames(vifTable)[ncol(vifTable)] <- "GVIF^(1/(2*Df))"
                    	printDataFrame(vifTable)
               	}
			cat("\n\n")
		}

		if (i == 1) {
			saveTable <- NULL
			tempDataX <- data.frame(rnamesX = as.numeric(rownames(data)), data)
			newName <- make.unique(c(names(data),paste(depVar[i],"pred", sep = "_")), sep = "")
			tempDataY <- data.frame(rnamesY = as.numeric(names(myfit[[i]]$fitted.values)), tempPred = myfit[[i]]$fitted.values)
			saveTable <- merge(tempDataX, tempDataY, by.x = "rnamesX", by.y = "rnamesY", all.x = TRUE)
			colnames(saveTable) <- c("rnamesX", newName)

		} else {
			newName <- make.unique(c(names(saveTable),paste(depVar[i],"pred", sep = "_")), sep = "")
			tempDataY <- data.frame(rnamesY = as.numeric(names(myfit[[i]]$fitted.values)), tempPred = myfit[[i]]$fitted.values)
			saveTable <- merge(saveTable, tempDataY, by.x = "rnamesX", by.y = "rnamesY", all.x = TRUE)
			colnames(saveTable) <- c(newName)
		}

		# SAVE RESIDUAL
		newName <- make.unique(c(names(saveTable),paste(depVar[i],"resid", sep = "_")), sep = "")
		tempDataY <- data.frame(rnamesY = as.numeric(names(myfit[[i]]$residuals)), tempResid = myfit[[i]]$residuals)
		saveTable <- merge(saveTable, tempDataY, by.x = "rnamesX", by.y = "rnamesY", all.x = TRUE)
		colnames(saveTable) <- newName

     	if (COOKS) {
          	newName <- make.unique(c(names(saveTable),paste(depVar[i],"cooks", sep = "_")), sep = "")
          	tempDataY <- data.frame(rnamesY = as.numeric(names(cooks.distance(myfit[[i]]))), tempCooks = cooks.distance(myfit[[i]]))
			saveTable <- merge(saveTable, tempDataY, by.x = "rnamesX", by.y = "rnamesY", all.x = TRUE)
		    	colnames(saveTable) <- newName
		}
		if (leverage) {
			tempDataY <- data.frame(rnamesY = as.numeric(names(hatvalues(myfit[[i]]))), tempHat = hatvalues(myfit[[i]]))
			newName <- make.unique(c(names(saveTable),paste(depVar[i],"leverage", sep = "_")), sep = "")
			saveTable <- merge(saveTable, tempDataY, by.x = "rnamesX", by.y = "rnamesY", all.x = TRUE)
			colnames(saveTable) <- newName
		}
	}
     
	if (is.null(saveTable)) { return(invisible(list(data = data, modelFit = myfit))) 
	} else { 
		saveTable <- saveTable[,-I(match("rnamesX", names(saveTable)))]
		return(invisible(list(data = saveTable, modelFit = myfit))) 
	}
}