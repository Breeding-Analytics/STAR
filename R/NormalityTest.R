# ----------------------------------------------------------------------------------------
# NormalityTest : Function for performing test for normality
# Created by: Alaine A. Gulles 02.24.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 08.15.2013 

# ----------------------------------------------------------------------------------------

NormalityTest <- function(data, var, grp = NULL, method = c("swilk")) UseMethod("NormalityTest")

NormalityTest.default <- function(data, var, grp = NULL, method = c("swilk")) {
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

	if (!is.character(var)) { stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
	if (!is.null(grp)) { 
		if (!is.character(grp)) { stop(paste("The object 'var' should be a character vector.", sep = "")) }
		if (any(is.na(match(grp, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
		for (i in (1:length(grp))) { if (!is.factor(tempData[,grp[i]])) { tempData[,grp[i]] <- factor(tempData[,grp[i]]) }}
	}

	procedure = c("swilk", "sfrancia", "ks", "cramer", "anderson")
	procedure.command <- c("shapiro.test", "sf.test", "lillie.test", "cvm.test", "ad.test")
	procedure.label <- c("Shapiro-Wilk", "Shapiro-Francia", "Kolmogorov-Smirnov", "Cramer-von Mises", "Anderson-Darling")
	method <- procedure[match(method, procedure)]
	if (all(is.na(method))) stop("The 'method' specified is invalid. Valid values for 'method' are 'swilk', 'sfrancia', 'ks', 'cramer', 'anderson'")

	method <- procedure[na.omit(match(method, procedure))]
	method.command <- procedure.command[match(method, procedure)]
	method.label <- procedure.label[match(method, procedure)]	

	normTable <- NULL

	if (is.null(grp)) {
		for (i in (1:length(method.command))) {
		     if (method[i] == "swilk") probSign <- "<" else probSign <- ">"
			for (j in (1:length(var))) {
				enough.sample <- TRUE
				if (method[i] == procedure[1])  { if (length(na.omit(tempData[,var[j]])) < 3 || length(na.omit(tempData[,var[j]])) > 5000)  enough.sample <- FALSE}
				if (method[i] == procedure[2])  { if (length(na.omit(tempData[,var[j]])) < 5 || length(na.omit(tempData[,var[j]])) > 5000)  enough.sample <- FALSE}
				if (method[i] == procedure[3])  { if (length(na.omit(tempData[,var[j]])) <= 5)  enough.sample <- FALSE}
				if (method[i] == procedure[4])  { if (length(na.omit(tempData[,var[j]])) <= 7)  enough.sample <- FALSE}
				if (method[i] == procedure[5])  { if (length(na.omit(tempData[,var[j]])) <= 7)  enough.sample <- FALSE}
				if (enough.sample) {
					result <- eval(parse(text = paste(method.command[i], "(tempData[,'",var[j],"'])", sep = "")))
					normTable <- rbind(normTable, data.frame(VARIABLE = var[j], PROCEDURE = method.label[i], 
                                                                  STATISTIC = names(result[[1]]), 
                                                                  VALUE = result[[1]][[1]], 
                                                                  PROB = paste("Pr(", probSign, " ",names(result[[1]]),")",sep = ""), 
                                                                  PVALUE = result[[2]]))
				}
			}
		}

          if (is.null(normTable)) { stop ("Number of observations is not enough to perform the test for normality.") }
          if (length(method) == 1) {
		     colnames(normTable) <- c("Variable", "Method", "Statistic", paste(levels(normTable[1,"STATISTIC"])[1], "Value"), "Prob", levels(normTable[1,"PROB"]))
		     normTable <- normTable[-I(c(3, 5))]
		} else { colnames(normTable) <- c("Variable", "Method", "Statistic", "Value","Prob", "p Value") }
          
	}

	if (!is.null(grp)) {
		for (i in (1:length(method.command))) {
               if (method[i] == "swilk") probSign <- "<" else probSign <- ">"
			for (j in (1:length(var))) {
				for (k in (1:length(grp))) {
					enough.sample <- TRUE
					tempData2 <- na.omit(data.frame(tempData[,var[j]], tempData[,grp[k]]))
					colnames(tempData2) <- c(var[j], grp[k])
					grpLevel <- levels(tempData[,grp[k]])
					if (method[i] == procedure[1])  { if (any(tapply(tempData2[,1], tempData2[,2], length) < 3) || any(tapply(tempData2[,1], tempData2[,2], length) > 5000))  { enough.sample <- FALSE }}
					if (method[i] == procedure[2])  { if (any(tapply(tempData2[,1], tempData2[,2], length) < 5) || any(tapply(tempData2[,1], tempData2[,2], length) > 5000))  { enough.sample <- FALSE }}
					if (method[i] == procedure[3])  { if (any(tapply(tempData2[,1], tempData2[,2], length) <= 5))  { enough.sample <- FALSE }}
					if (method[i] == procedure[4])  { if (any(tapply(tempData2[,1], tempData2[,2], length) <= 7))  { enough.sample <- FALSE }}
					if (method[i] == procedure[5])  { if (any(tapply(tempData2[,1], tempData2[,2], length) <= 7))  { enough.sample <- FALSE }}
					if (enough.sample) {
						result <- tapply(tempData2[,1], tempData2[,2], method.command[i])
						for (l in (1:length(result))) { 
							normTable <- rbind(normTable, data.frame( GRP = grp[k], 
														LEVEL = grpLevel[l], 
														VARIABLE = var[j], 
														PROCEDURE = method.label[i], 
                                                                      STATISTIC = names(result[[l]][[1]]),       
														VALUE = result[[l]][[1]][[1]], 
                                                                      PROB = paste("Pr(", probSign, " ",names(result[[l]][[1]]),")",sep = ""),       
														PVALUE = result[[l]][[2]])) 
						}
					}
				}
			}
		}
          
		if (is.null(normTable)) { stop ("Number of observations is not enough to perform the test for normality.") }
		
		normTable <- normTable[order(normTable$GRP, normTable$VARIABLE, normTable$LEVEL),]
		rownames(normTable) <- c(1:nrow(normTable))
		
          if (length(method) == 1) {
               colnames(normTable) <- c("Grp", "Level", "Variable", "Method", "Statistic", paste(levels(normTable[1,"STATISTIC"])[1], "Value"), "Prob", levels(normTable[1,"PROB"]))
               normTable <- normTable[-I(c(5, 7))]
          } else { colnames(normTable) <- c("Grp", "Level", "Variable", "Method", "Statistic", "Value", "Prob", "p Value") }
		
	}

	if (!is.null(normTable)) { 
		if (length(var) == 1) { tempVar <- c(rep(var, nrow(normTable))) }
		options(width = 5000)
          if (length(method) == 1) {  cat("Test for Normality","\n")     
          } else {  cat("Tests for Normality","\n") }
		printDataFrame(normTable)
	}

	return(invisible(normTable))
} ### end stmt -- NormalityTest 