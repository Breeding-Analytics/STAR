# -----------------------------------------------------------------
# BALANCED INCOMPLETE BLOCK DESIGN
# File Created by: Alaine A. Gulles 04.13.2011
#                  for International Rice Research Institute
# File Modified by: Alaine A. Gulles 10.14.2013
# Reference: agricolae package
# -----------------------------------------------------------------

BIBDTest <- function(data, respvar, trmt, block, set = NULL, method = NULL, descriptive = FALSE, normality = FALSE, homogeneity = FALSE, alpha = 0.05, outputPath = NULL) UseMethod("BIBDTest")

BIBDTest.default <- function(data, respvar, trmt, block, set = NULL, method = NULL, descriptive = FALSE, normality = FALSE, homogeneity = FALSE, alpha = 0.05, outputPath = NULL) {

     options(width = 5000, show.signif.star = FALSE)
     if (is.character(data)) { data <- eval(parse(text = data)) } 
     
	# define trmt as factors
	if (is.numeric(data[,trmt])) { data[,trmt] <- factor(data[,trmt]) 
	} else {
	     if (suppressWarnings(all(!is.na(as.numeric(as.character(factor(trimStrings(data[,trmt])))))))) {
	          data[,trmt] <- factor(as.numeric(data[,trmt]))
	     } else { data[,trmt] <- factor(trimStrings(data[,trmt])) }   
	}
     
     # define block as factor
	if (is.numeric(data[,block])) { data[,block] <- factor(data[,block]) 
	} else {
	     if (suppressWarnings(all(!is.na(as.numeric(as.character(factor(trimStrings(data[,block])))))))) {
	          data[,block] <- factor(as.numeric(data[,block]))
	     } else { data[,block] <- factor(trimStrings(data[,block])) }   
	}
     
	#procedure <- c("lsd", "tukey", "snk", "duncan")
	#method <- procedure[na.omit(match(method, procedure))]
     #if (length(method) == 0) method <- NULL

	LSDTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, trmtlevels) {
		Tprob <- qt(1 - alpha/2, dfError)
		thestring <- paste(paste(capture.output(print(output1 <- order.group(trmtlevels,
		      	      TrmtMean.adj, rep(r, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), StdErr.diff))), collapse = "\n"), collapse = "\n\n")
		output1[,"std.err"] <- NULL
		colnames(output1)[1] <- trmt
          output1 <- output1[order(output1[,1]),]
		return(list(method = "Least Significant Difference (LSD) Test", cval = Tprob * StdErr.diff, stderr = StdErr.diff, pw = output1, alpha = alpha))				
	}

	TukeyTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, trmtlevels) {
		Tprob <- qtukey(1 - alpha, a, dfError)
		thestring <- paste(paste(capture.output(print(output1 <- order.group(trmtlevels,
				      TrmtMean.adj, rep(r, a), MSError, Tprob, StdErr.adjtrtmean,
					(k /(lambda * a))/2, StdErr.diff))), collapse = "\n"), collapse = "\n\n")
		output1[,"std.err"] <- NULL
		colnames(output1)[1] <- trmt
		output1 <- output1[order(output1[,1]),]
		return(list(method = "Tukey's Honestly Significant Difference (HSD) Test", 
				cval = Tprob * StdErr.diff/sqrt(2), stderr = StdErr.diff/sqrt(2), pw = output1, alpha = alpha))
	}

	DuncanTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, trmtlevels){
		Tprob <- qtukey((1 - alpha)^(1:(a-1)), 2:a, dfError)
		duncan <- Tprob * StdErr.diff/sqrt(2)
		critical.range <- rbind(tabular.value = Tprob, critical.value = duncan)
          colnames(critical.range) <- 2:a
		thestring <- paste(paste(capture.output(print(output1 <- order.group(trmtlevels,
			            TrmtMean.adj, rep(r, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), 2, dfError, alpha, StdErr.diff/sqrt(2)))), collapse = "\n"), collapse = "\n\n")
		output1[,"std.err"] <- NULL
		colnames(output1)[1] <- trmt
		output1 <- output1[order(output1[,1]),]
		return(list(method = "Duncan's Multiple Range test", crange = critical.range, stderr = StdErr.diff/sqrt(2), pw = output1, alpha = alpha))
	}

	SNKTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, trmtlevels){
		Tprob <- qtukey((1 - alpha), 2:a, dfError)
		SNK <- Tprob * StdErr.diff/sqrt(2)
		critical.range <- rbind(tabular.value = Tprob, critical.value = SNK)
		colnames(critical.range) <- 2:a
		thestring <- paste(paste(capture.output(print(output1 <- order.group(trmtlevels,
			           TrmtMean.adj, rep(r, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), 2, dfError, alpha, StdErr.diff/sqrt(2)))), collapse = "\n"), collapse = "\n\n")
		output1[,"std.err"] <- NULL
		colnames(output1)[1] <- trmt
		output1 <- output1[order(output1[,1]),]
		return(list(method = "Student Newman Keuls (SNK) Test", crange = critical.range, stderr = StdErr.diff/sqrt(2), pw = output1, alpha = alpha))
	}

	anova.table <- list()
	trmtStat <- NULL
	blkStat <- NULL
	tempStat <- NULL
     rvWithSigEffect <- NULL
     
     if (is.null(set)) {
          set <- make.unique(c(names(data), "setVar"))[length(make.unique(c(names(data), "setVar")))]
          data <- cbind(data, setVar = "1")
          addSet <- TRUE
     } else { 
          addSet <- FALSE 
          data[,set] <- factor(data[,set])
     }
     tempAllData <- data
     
	cat("Analysis of Variance","\n", "Balanced Incomplete Block" ,"\n\n", sep = "")

	for (i in (1:length(respvar))) {
	     width1 <- 37 + nchar(respvar[i])
	     cat(paste(rep("=", width1), collapse = ""), "\n")
	     cat("ANALYSIS FOR RESPONSE VARIABLE:", respvar[i],"\n")
	     cat(paste(rep("=", width1), collapse = ""), "\n\n")

          anova.table[[i]] <- list()
	     anova.table[[i]]$respvar <- respvar[i]
	     anova.table[[i]]$set <- list()
	     tempNewData <- NULL

	     for (z in (1:nlevels(data[,match(set, names(data))]))) {
	          estimateData <- FALSE
	          setLabel <- levels(data[,match(set, names(data))])[z]
	          anova.table[[i]]$set[[z]] <- list()
	          anova.table[[i]]$set[[z]]$setValue <- setLabel
               
	          if (!addSet) {
	               width2 <- 5 + nchar(set) + nchar(setLabel)
	               cat(paste(rep("-", width2), collapse = ""), "\n")
	               cat(set, "=", setLabel, "\n")
	               cat(paste(rep("-", width2), collapse = ""), "\n\n")
	               rm(width2)
	          } ### end stmt -- if (!addSet)  
               
               # subset the data set by set variable
	          origData <- tempAllData[tempAllData[,set] == levels(tempAllData[,match(set, names(tempAllData))])[z],]
               if (length(na.omit(origData[,respvar[i]])) == 0) {
                    cat("Error: The data set does not contain any observations.\n\n")
                    next
               }
               
               # convert the block and trmt as factor in the subset data
               origData[,block] <- factor(origData[,block])
	          origData[,trmt] <- factor(origData[,trmt])
                              
	          tempData <- origData
               
               # compute the levels of the treatment and the block, rep, k and lambda
	          a <- nlevels(tempData[,trmt]) # -- trmt levels
	          b <- nlevels(tempData[,block]) # -- blk levels
	          r <- unique(table(tempData[,trmt]))
	          k <- unique(table(tempData[,block]))
	          lambda <- (r * (k - 1))/(a - 1)
               
               # determine if the data set is balanced
	          capture.output(numObsTrmt <- DescriptiveStatistics(data = tempData, var = respvar[i], grp = trmt)[,"N_NonMissObs"])
               if (!all(numObsTrmt == r)) {
                    cat("Error: Cannot perform ANOVA for response variable '", respvar[i], "'. The data set is unbalanced.",sep = "", "\n\n")
                    next
               }
               
	          # determine if the data set is balanced
	          capture.output(numObsBlk <- DescriptiveStatistics(data = tempData, var = respvar[i], grp = block)[,"N_NonMissObs"])
	          if (!all(numObsBlk == k)) {
	               cat("Error: Cannot perform ANOVA for response variable '", respvar[i], "'. The data set is unbalanced.",sep = "", "\n\n")
	               next
	          }
               
               # perform analysis of variance for adjusted treatment means
	          withError <- try(result1 <- lm(formula(paste(respvar[i], "~", block, "+", trmt)), data = tempData), silent = TRUE)
               
               if (class(withError) == "try-error") {
                    msg <- trimStrings(strsplit(withError, ":")[[1]])
                    msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
                    msg <- gsub("\"", "", msg)
                    cat("Error: ", msg, sep = "", "\n\n")
                    next
               } 
                              
	          tempAnova1 <- anova(result1)
	          rownames(tempAnova1)[1] <- paste(rownames(tempAnova1)[1], "(unadj)", sep = "")
	          rownames(tempAnova1)[2] <- paste(rownames(tempAnova1)[2], "(adj)", sep = "")
	          rownames(tempAnova1)[nrow(tempAnova1)] <- "Error"
	          tempAnova1[nrow(tempAnova1)+1, 1:2] <- c(sum(tempAnova1[,1]), sum(tempAnova1[,2]))
	          rownames(tempAnova1)[nrow(tempAnova1)] <- "Total"
	          tempAnova1[1,3:5] <- NA
               
	          # perform analysis of variance for adjusted block means
	          result2 <- lm(formula(paste(respvar[i], "~", trmt, "+", block)), data = tempData)
	          tempAnova2 <- anova(result2)     
	          rownames(tempAnova2)[1] <- paste(rownames(tempAnova2)[1], "(unadj)", sep = "")
	          rownames(tempAnova2)[2] <- paste(rownames(tempAnova2)[2], "(adj)", sep = "")
	          tempAnova2[1,3:5] <- NA
               
	          anova.table[[i]]$set[[z]]$ANOVA.adjtrmt <- suppressWarnings(anova(result1))
	          anova.table[[i]]$set[[z]]$ANOVA.adjblk <- suppressWarnings(anova(result2))
               aovTable <- rbind(tempAnova1[2,], tempAnova2[1,], tempAnova1[1,], tempAnova2[2,], tempAnova1[3,], tempAnova1[4,])
	          anova.table[[i]]$set[[z]]$ANOVA.table <- aovTable
               #anova.table[[i]] <- rbind(tempAnova1[2,], tempAnova2[1,], tempAnova1[1,], tempAnova2[2,], tempAnova1[3,], tempAnova1[4,])
	          dfError <- tempAnova1[nrow(tempAnova1)-1,1]
	          MSError <- tempAnova1[nrow(tempAnova1)-1,3]
	          anova.table[[i]]$set[[z]]$dfError <- dfError
	          anova.table[[i]]$set[[z]]$MSError <- MSError
               
	          # display summary information
	          ClassInformation(tempData[,c(block, trmt, respvar[i])], respvar = respvar[i])
	          cat("\n\n")
               
               # display descriptive statistics
	          if (descriptive) {
	               DescriptiveStatistics(data, var = respvar[i], statistics = c("n", "mean", "sd", "min", "max"))
	               cat("\n")
	          }
               
               # create the residual data and predicted values
               residNfittedData <- NULL
	          residNfittedData <- data.frame(PredictedValues(result1), resid(result1))
               newNames <- make.unique(c(names(tempData), paste(respvar[i],"_pred", sep = ""), paste(respvar[i],"_resid", sep = "")), sep = "")
               colnames(residNfittedData) <- newNames[(length(newNames)-1):length(newNames)]
               
               
	          # -- CREATE THE DIAGNOSTIC PLOT -- #
	          if (!is.null(outputPath)) {
	               if (addSet) { png(filename = paste(outputPath, "BIBDDiagPlot_", respvar[i], ".png", sep = ""), width = 480, height = 240)
	               } else { png(filename = paste(outputPath, "BIBDDiagPlot_", setLabel,"_", respvar[i], ".png", sep = ""), width = 480, height = 240) }
	               params <- par(mfrow = c(1,2), bg = "white")
	               if (addSet) { plot(residNfittedData[,(ncol(residNfittedData)-1)], residNfittedData[,ncol(residNfittedData)], main = paste("Residual vs Fitted:\n", respvar[i], sep = ""), xlab = "Fitted Values", ylab = "Residuals") 
	               } else { plot(residNfittedData[,(ncol(residNfittedData)-1)], residNfittedData[,ncol(residNfittedData)], main = paste("Residual vs Fitted:\n", respvar[i], " for ", set, " = ", setLabel, sep = ""), xlab = "Fitted Values", ylab = "Residuals") }
	               qqnorm(residNfittedData[,ncol(residNfittedData)])
	               qqline(residNfittedData[,ncol(residNfittedData)])
	               par(params)
	               # params <- par(mfrow = c(2,2), bg = "white")
	               # plot(result)
	               # par(params)
	               dev.off()
	          }
               
               # perfomr normality test and test for homogeneity of variances
               
	          if (homogeneity || normality) {
	               tempData2 <- tempData[!is.na(tempData[,respvar[i]]),]
	               assumpData <- data.frame(tempData2, residNfittedData[ncol(residNfittedData)])
	               
	               # --- PRINTING RESULTS OF TEST FOR HOMOGENEITY OF VARIANCES --- #
	               if (homogeneity) {
	                    capture.output(bartlett.result <- HeteroskedasticityTest(data = assumpData, var = names(assumpData)[ncol(assumpData)], grp = trmt, method = c("bartlett")))
	                    cat("Test for Homogeneity of Variances\n")
	                    printDataFrame(bartlett.result[,3:ncol(bartlett.result)])
	                    cat("\n")
	                    rm(bartlett.result)
	               }
	               
	               # --- PRINTING RESULT OF SHAPIRO WILK TEST --- #
	               if (normality) {
	                    if (nrow(assumpData) >= 3 && nrow(assumpData) <= 5000) {
	                         NormalityTest(data = assumpData, var = names(assumpData)[ncol(assumpData)], grp = NULL, method = c("swilk"))
	                         cat("\n")
	                    }
	               }
	               rm(assumpData)
	          }
	          
               # print the anova table
	          cat("ANOVA TABLE\nResponse Variable: ", respvar[i], "\n", sep = "")
	          printAOVTable(aovTable)
               cat("\n")
	          
	          # computation of summary statistics for adjusted and unadjusted trmt means
	          blk.sum <- tapply(tempData[,respvar[i]], tempData[,block], sum, na.rm = TRUE)
	          trt.sum <- tapply(tempData[,respvar[i]], tempData[,trmt], sum, na.rm = TRUE)
	          trt.mean <- tapply(tempData[,respvar[i]], tempData[,trmt], mean, na.rm = TRUE)
	          cmb.mean <- tapply(tempData[,respvar[i]], list(tempData[,trmt], tempData[,block]), mean, na.rm = TRUE)
	          num <- (!is.na(cmb.mean)) %*% matrix(blk.sum)
	          trmtTotal.adj <- matrix(trt.sum) - (as.numeric(num)/k)
	          SSTrmt.adj <- (colSums(trmtTotal.adj * trmtTotal.adj) * k)/(lambda * a)
	          MSTrmt.adj <- SSTrmt.adj/(a - 1)
	          StdErr.diff <- sqrt((2 * k * MSError)/(lambda * a))
	          StdErr.adjtrteff <- sqrt((k * MSError)/(lambda * a))
	          StdErr.adjtrtmean <- sqrt((MSError * (1 + ((k * r * (a - 1))/(lambda * a))))/(a * r))
	          #sqrt((MSError * (1 + ((k * numObsTrmt * (a - 1))/(lambda * a))))/(a * numObsTrmt))
	          TrmtEffect.adj <- trmtTotal.adj * (k /(lambda * a))    # --- TrmtEffect.adj = intrablock estimate
	          interblk.est <- (num - (k * r * mean(tempData[,respvar[i]], na.rm = TRUE)))/(r - lambda)
	          TrmtMean.adj <- mean(tempData[,respvar[i]], na.rm = TRUE) + TrmtEffect.adj
	          
               # print summary statistics
	          cat("Table of Means\n", sep = "")
	          tempTable <- data.frame(RespVar = respvar[i], TrmtMean = trt.mean, AdjTrmtMean = TrmtMean.adj, AdjTrmtEffect = TrmtEffect.adj)
	          tempTable <- data.frame(rownames(tempTable), tempTable)
	          colnames(tempTable)[1] <- trmt
	          printDataFrame(tempTable)
	          trmtStat <- rbind(trmtStat, cbind(setLabel, tempTable))
                              
	          cat("\nStandard Error\n", sep = "")
	          tempTable <- data.frame(Diff = StdErr.diff, AdjTrmtEffect = StdErr.adjtrteff, AdjTrmtMean = StdErr.adjtrtmean)
	          printDataFrame(tempTable)
	          tempStat <- rbind(tempStat, cbind(setLabel, data.frame(RespVar = respvar[i], tempTable)))
	                         
               cat("\nInterblock Estimate\n", sep = "")
	          tempTable <- data.frame(RespVar = respvar[i], t(interblk.est))
               colnames(tempTable)[2:ncol(tempTable)] <- levels(tempData[,trmt]) 
	          printDataFrame(tempTable)
	          tempTable <- data.frame(setLabel, RespVar = respvar[i], levels(tempData[,trmt]), interblk.est)
               names(tempTable)[3] <- trmt
	          blkStat <- rbind(blkStat, tempTable)
               cat("\n")
               
               # determine if the adjusted treatment means is significant
	          if (aovTable[match(paste(trmt, "(adj)", sep = ""), rownames(aovTable)),5] < alpha) {
                    rvWithSigEffect <- c(rvWithSigEffect, respvar[i])
                    anova.table[[i]]$set[[z]]$sigTrmt <- TRUE
                    anova.table[[i]]$set[[z]]$pw <- list()
                    anova.table[[i]]$set[[z]]$pw$alpha <- alpha
                    anova.table[[i]]$set[[z]]$pw$LSD <- LSDTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, levels(tempData[,trmt])) 
                    anova.table[[i]]$set[[z]]$pw$HSD <- TukeyTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, levels(tempData[,trmt])) 
                    anova.table[[i]]$set[[z]]$pw$SNK <- SNKTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, levels(tempData[,trmt]))
                    anova.table[[i]]$set[[z]]$pw$duncan <- DuncanTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff, alpha, a, lambda, r, k, levels(tempData[,trmt])) 
                    
                    if (is.null(method)) {
                         if (nlevels(data[,trmt]) <= 5) { printBIBDPairwise(anova.table[[i]]$set[[z]]$pw, dfError, MSError, method = "LSD", trmtLabel = trmt)
                         } else { printBIBDPairwise(anova.table[[i]]$set[[z]]$pw, dfError, MSError, method = "HSD", trmtLabel = trmt) }
                    } else {
                         printBIBDPairwise(anova.table[[i]]$set[[z]]$pw, dfError, MSError, method, trmtLabel = trmt)
                    }
	          } else { 
	               anova.table[[i]]$set[[z]]$sigTrmt <- FALSE
	               anova.table[[i]]$set[[z]]$pw <- NULL
	          }
               
               # combining the data set
               if (is.null(tempNewData)) { tempNewData <- cbind(tempData, residNfittedData) 
               } else { 
                    if (ncol(cbind(tempData, residNfittedData)) > ncol(tempNewData)) {
                         tempName <- names(cbind(tempData, residNfittedData))[-I(match(names(tempNewData), names(cbind(tempData, residNfittedData))))]
                         tempNewData <- cbind(tempNewData, rep("obs", nrow(tempNewData)))
                         names(tempNewData)[ncol(tempNewData)] <- tempName
                    } else {
                         if (ncol(cbind(tempData, residNfittedData)) < ncol(tempNewData)) {
                              tempName <- names(tempNewData)[-I(match(names(cbind(tempData, residNfittedData)), names(tempNewData)))]
                              tempData <- cbind(tempData, rep("obs", nrow(tempData)))
                              names(tempData)[ncol(tempData)] <- tempName
                         }     
                    }
                    tempNewData <- rbind(tempNewData, cbind(tempData, residNfittedData))
               }
          } ### end stmt -- for (z in (1:nlevels(data[,match(set, names(data))])))
          
	     if (!is.null(tempNewData)) tempAllData <- merge(tempAllData, tempNewData, all = TRUE)
     
	} ### end stmt -- for (i in (1:length(respvar))) 

     rownames(tempStat) <- 1:nrow(tempStat)
	rownames(trmtStat) <- 1:nrow(trmtStat)
	rownames(blkStat) <- 1:nrow(blkStat)
     if (!is.null(rvWithSigEffect)) { rvWithSigEffect <- unique(rvWithSigEffect) }
     if (addSet) {
          tempAllData[,"setVar"] <- NULL
          tempStat[,"setLabel"] <- NULL
          trmtStat[,"setLabel"] <- NULL
          blkStat[,"setLabel"] <- NULL
          set <- NULL
     } else {
          names(tempStat)[match("setLabel", names(tempStat))] <- set
          names(trmtStat)[match("setLabel", names(trmtStat))] <- set
          names(blkStat)[match("setLabel", names(blkStat))] <- set
     }
     
     return(invisible(list(data = tempAllData, result = anova.table, trmtStatistic = trmtStat, stdErr = tempStat, InterBlkEst = blkStat, rvWithSigEffect = rvWithSigEffect, set = set)))
}
