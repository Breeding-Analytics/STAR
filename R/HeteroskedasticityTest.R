# -------------------------------------------------------------------------------------
# RCropStat Beta Version: Function for DATA MENU
# -------------------------------------------------------------------------------------
# HeteroskedasticityTest: Perform test for homogeneity of variances
# Created by: Alaine A. Gulles 03.16.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 02.16.2012 
# -------------------------------------------------------------------------------------
# ARGUMENT: var   = numeric data frame or a character vector containing the variable 
#			  names needed for the analysis
#           grp = a numeric or character data frame or a character vector containing
#			the grouping variable names needed for the analysis
#		data = data frame containing the variables needed in the analysis
#		method = a character string specifying the method to be used, must at least
#			   one of the followng "bartlett" (default) and/or "levene"
# -------------------------------------------------------------------------------------

HeteroskedasticityTest <- function(data, var, grp, method = c("bartlett")) UseMethod("HeteroskedasticityTest")

HeteroskedasticityTest.default <- function(data, var, grp, method = c("bartlett")) {
     if (is.character(data)) { 
          nameData <- data
          if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
          tempData <- eval(parse(text = data))
     } else {
          if (is.data.frame(data)) {
               nameData <- paste(deparse(substitute(data)))	
               tempData <- data
          } else {
               stop ("The argument should either be a data frame or a character string indicating the name of the data frame.")
          }
     }
     if (is.null(var) || is.null(grp)) { stop("The object 'var' and 'grp' should be a non-null character vector.") }
	if (!is.character(var)) { stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
	if (!is.character(grp)) { stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(grp, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }

	procedure = c("bartlett", "levene")
	procedure.command = c("bartlett.test", "leveneTest")
	method <- procedure[match(method, procedure)]
	if (all(is.na(method))) stop("The 'method' specified is invalid. Valid values for 'method' are 'bartlett' and 'levene'")
	method <- procedure[na.omit(match(method, procedure))]
	method.command <- procedure.command[match(method, procedure)]

	hovTable <- NULL

	for (i in (1:length(method))) {
		for (j in (1:length(var))) {
			for (k in (1:length(grp))) {
                    if (!is.factor(tempData[,grp[k]])) { tempData[,grp[k]] <- factor(tempData[,grp[k]]) }
				command <- paste(method.command[i], "(tempData[,'", var[j],"'], tempData[,'", grp[k],"'])", sep = "")
				result <- eval(parse(text = command))
				if (method[i] == "levene") hovTable <- rbind(hovTable, data.frame(GRP = grp[k], VARIABLE = var[j],
                                                                                      PROCEDURE = paste(toupper(substring(method[i],1,1)),substring(method[i],2), sep = ""), 
                                                                                      DF = result[1,1], STATISTIC = "F", VALUE = result[1,2], 
                                                                                      PROB = colnames(result)[3], PVALUE = result[1,3]))
				if (method[i] == "bartlett") hovTable <- rbind(hovTable, data.frame(GRP = grp[k], VARIABLE = var[j], 
                                                                PROCEDURE = paste(toupper(substring(method[i],1,1)),substring(method[i],2), sep = ""), 
                                                                DF = result[[2]], STATISTIC = "Chisq",
                                                                VALUE = result[[1]], PROB = "Pr(>Chisq)", PVALUE = result[[3]]))
			
			}
		}
	}
     
     hovTable <- hovTable[order(hovTable$GRP, hovTable$VARIABLE),]

	if (!is.null(hovTable)) { 
		rownames(hovTable) <- c(1:nrow(hovTable))
          if (length(method) == 1) {
               colnames(hovTable) <- c("Grp", "Variable", "Method", "DF", "Statistic", paste(levels(hovTable[1,"STATISTIC"])[1], "Value"), "Prob", levels(hovTable[1,"PROB"])[1])
               hovTable <- hovTable[-I(c(5,7))]
          } else {
               colnames(hovTable) <- c("Grp", "Variable", "Method", "DF", "Statistic", "Value", "Prob", "p Value")
          }
		options(width = 5000)
          if (length(method) == 1) {  cat("Test for Homogeneity of Variances","\n")
          } else {  cat("Tests for Homogeneity of Variances","\n")  }
		
		#print(hovTable, right = FALSE, row.names = FALSE)
		printDataFrame(hovTable)
	}
	return(invisible(hovTable))
} ### end stmt -- HeteroskedasticityTest 