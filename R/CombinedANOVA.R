# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# combAOVTest: Functions for performing ANOVA
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.05.2013
# Note: Include remarks when data is unbalanced, when no balanced data can be generated
# -------------------------------------------------------------------------------

combAOVTest <- function(data, design, respvar, factor1, factor2 = NULL, factor3 = NULL, 
                        factor4 = NULL, rep1 = NULL, rep2 = NULL, set, descriptive = FALSE, 
                        normality = FALSE, homogeneity = FALSE, pwTest = NULL, pwVar = NULL, 
                        contrastOption = NULL, sig = 0.05, outputPath = NULL) UseMethod("combAOVTest")


combAOVTest.default <- function(data, design, respvar, factor1, factor2 = NULL, factor3 = NULL, 
                        factor4 = NULL, rep1 = NULL, rep2 = NULL, set, descriptive = FALSE, 
                        normality = FALSE, homogeneity = FALSE, pwTest = NULL, pwVar = NULL, 
                        contrastOption = NULL, sig = 0.05, outputPath = NULL) {

     if (is.character(data)) { 
          nameData <- data
          data <- eval(parse(text = data))
     } else { nameData <- paste(deparse(substitute(data))) }
     
     if (!is.data.frame(data)) { stop("The object 'data' should be a data frame.") }
     
     availableDesign <- c("CRD", "RCBD", "LSD", "SplitCRD", "SplitRCBD", "SplitLSD", "Strip", "Split2CRD", "Split2RCBD", "Split2LSD", "Strip-Split",     "Split3CRD", "Split3RCBD", "Split3LSD", "Strip-Split2")
     if(is.na(match(design, availableDesign))) {
          stop("Design must be one of the following:\n'CRD', 'RCBD', 'LSD',\n'SplitCRD', 'SplitRCBD', 'SplitLSD','Strip',\n'Split2CRD', 'Split2RCBD', 'Split2LSD', 'Strip-Split',\n'Split3CRD', 'Split3RCBD', 'Split3LSD', 'Strip-Split2'")
     }
     
     designChoice <- match(design, availableDesign) ## design code
     designTitle <- c("Completely Randomized Design",
                      "Randomized Complete Block Design",
                      "Latin Square Design", 
                      "Split Plot in Completely Randomized Design", 
                      "Split Plot in Randomized Complete Block Design", 
                      "Split Plot in Latin Square Design", 
                      "Strip Plot Design", 
                      "Split-Split Plot in Completely Randomized Design", 
                      "Split-Split Plot in Randomized Complete Block Design", 
                      "Split-Split Plot in Latin Square Design",  
                      "Strip-Split Plot Design", 
                      "Split-Split-Split Plot in Completely Randomized Design", 
                      "Split-Split-Split Plot in Randomized Complete Block Design", 
                      "Split-Split-Split Plot in Latin Square Design", 
                      "Strip-Split-Split Plot Design")
     
     allFactor <- c(factor1, factor2, factor3, factor4)
     prev.option <- options()$show.signif.stars
     options(show.signif.stars = FALSE, width = 5000)
     options(expressions = 500000)
     tempAnova <- list()
     residNfittedData <- NULL
     pwOption <- list()
     aovresult <- list()
     rvWithSigEffect <- NULL
     
     # define as factors
     for (i in (1:length(allFactor))) { data[,allFactor[i]]     <- factor(data[,allFactor[i]])	 }
     if (!is.null(rep1)) { data[,rep1] <- factor(data[,rep1]) }
     if (!is.null(rep2)) { data[,rep2] <- factor(data[,rep2]) }
     data[,set] <- factor(data[,set])
     
     if (nlevels(data[,set]) < 2) { stop(paste("The variable '",set,"' should have at least two levels to perform the analysis.", sep = "")) }
     
     tempData <- data
     #tempNewData <- NULL
     cat("Combined Analysis of Variance","\n",designTitle[designChoice],"\n\n", sep = "")
     
     for (i in (1:length(respvar))) {
          width1 <- 37 + nchar(respvar[i])
          cat(paste(rep("=", width1), collapse = ""), "\n")
          cat("ANALYSIS FOR RESPONSE VARIABLE:", respvar[i],"\n")
          cat(paste(rep("=", width1), collapse = ""), "\n\n")
          
          origData <- tempData
          #origData <- data
          tempNewData <- NULL
          formula <- paste(respvar[i], " ~ ", paste(c(set, rep1, rep2, allFactor), collapse = ":", sep = ""), sep = "")
          
          # determine if balanced data
          if (designChoice != 1) {
               rawData <- tempData[,sort(match(c(respvar[i], allFactor, rep1, rep2, set), names(tempData)))]
               if (nrow(na.omit(rawData)) == 0) {
                    cat("Error: Data set has no observation.\n\n", sep = "")
                    next
               }
               if (designChoice != 3 && designChoice != 6 && designChoice != 10 && designChoice != 14) {
                    if (is.list(replications(formula, rawData))) {
                         if (max(replications(formula, rawData)[[1]]) != 1) {
                              cat("ERROR: Cannot perform ANOVA for balanced data. The data set for response variable '", respvar[i],"' is unbalanced. \n\n", sep = "")
                              next
                         } 
                         for (j in (1:nlevels(tempData[,set]))) {
                              tempGenData <- GenerateBalanceData(tempData[tempData[, set] == levels(tempData[,set])[j],], respvar[i], allFactor, c(rep1, rep2) , design)
                              tempGenData[,set] <- levels(tempData[,set])[j]
                              tempNewData <- rbind(tempNewData, tempGenData)
                              #tempNewData <- rbind(tempNewData, GenerateBalanceData(tempData[tempData[, set] == levels(tempData[,set])[j],], respvar[i], allFactor, c(rep1, rep2) , design))
                         }
                         tempNewData[,set] <- factor(tempNewData[,set])
                         if ((1 - (length(na.omit(tempNewData[,respvar[i]]))/nrow(tempNewData))) > 0.10) {
                              cat("ERROR: Cannot perform ANOVA for balanced data for response variable '", respvar[i], "'. Too many missing values.\n\n", sep = "")
                              next
                         }
                         tempEstData <- NULL
                         for (j in (1:nlevels(tempData[,set]))) {
                              #tempEstData <- rbind(tempEstData, estMissData(design, data = tempData[tempData[, set] == levels(tempData[,set])[j],], respvar[i], factor1, factor2, factor3, factor4, rep1, rep2))
                              tempEstData <- rbind(tempEstData, estMissData(design, data = tempNewData[tempNewData[, set] == levels(tempNewData[,set])[j],], respvar[i], factor1, factor2, factor3, factor4, rep1, rep2))
                         }
                         tempData <- tempEstData
                         estimatedData <- TRUE
                    } else { estimatedData <- FALSE } 
               } else {
                    formula1 <- paste(respvar[i]," ~ ", paste(c(rep1, allFactor), collapse = ":", sep = ""), sep = "") 
                    if (is.list(replications(formula1, rawData))) {
                         if (max(replications(formula1, rawData)[[1]]) != 1) {
                              cat("ERROR: Cannot perform ANOVA for balance data. The data set for response variable '", respvar[i],"' is unbalanced.\n\n", sep = "")
                              next
                         } 
                         for (j in (1:nlevels(tempData[,set]))) {
                              tempNewData <- rbind(tempNewData, GenerateBalanceData(tempData[tempData[, set] == levels(tempData[,set])[j],], respvar[i], allFactor, c(rep1, rep2) , design))
                         }
                         if ((1 - (length(na.omit(tempNewData[,respvar[i]]))/nrow(tempNewData))) > 0.10) {
                              cat("ERROR: Cannot perform ANOVA for balanced data for response variable '", respvar[i], "'. Too many missing values.\n\n", sep = "")
                              next
                         }
                         tempEstData <- NULL
                         for (j in (1:nlevels(tempData[,set]))) {
                              tempEstData <- rbind(tempEstData, estMissData(design, data = tempData[tempData[, set] == levels(tempData[,set])[j],], respvar[i], factor1, factor2, factor3, factor4, rep1, rep2))
                         }
                         tempData <- tempEstData
                         estimatedData <- TRUE
                    } else { estimatedData <- FALSE }
               }
          } else { estimatedData <- FALSE }
          

          # -- PRINTING CLASS LEVEL INFORMATION -- #
          #ClassInformation(tempData[, c(set, rep1, rep2, factor1, factor2, factor3, factor4, respvar[i])], respvar = respvar[i])
          ClassInformation(data[, c(set, rep1, rep2, factor1, factor2, factor3, factor4, respvar[i])], respvar = respvar[i])
          cat("\n\n")
          
          # --- CHECKING OF ASSUMPTIONS --- #
          #if (homogeneity) {
          #     prelimCombAOV(model = buildAOVModel(design, respvar[i], factor1, factor2, factor3, factor4, rep1, rep2, set = NULL)[[1]], data =  tempData, set)
          #     cat("\n")
          #}
          
          
          # --- PRINTING DESCRIPTIVE STATISTICS --- #     
          if (descriptive) { 
               DescriptiveStatistics(data = tempData, var = respvar[i], grp = set, statistics = c("n", "mean", "sd", "min", "max"))
               if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
               cat("\n")
          }
          
          # -- BUILD THE MODEL -- #
          tempModel <- buildAOVModel(design, respvar[i], factor1, factor2, factor3, factor4, rep1, rep2, set)
          mymodel <- tempModel[[1]]
          mymodel2 <- tempModel[[2]]
          modelRHS <- tempModel[[3]]
          modelRHS2 <- tempModel[[4]]
          
          # --- PERFORM ANALYSIS OF VARIANCE --- #
          if (estimatedData) { tempresult <- summary(suppressWarnings(aov(formula(mymodel), data = origData)))     }
          result <- suppressWarnings(aov(formula(mymodel), tempData))
          aovresult[[i]] <- result
          
          if (estimatedData) {
               if (attr(summary(result), "class")[[1]] == "summary.aovlist") {
                    tempAnova[[i]] <- summary(result) 
                    numRow <- nrow(tempAnova[[i]][[length(tempAnova[[i]])]][[1]])
                    dfError <- tempresult[[length(tempresult)]][[1]][nrow(tempresult[[length(tempresult)]][[1]]),"Df"]
                    tempAnova[[i]][[length(tempAnova[[i]])]][[1]][numRow,"Df"] <- dfError
                    tempAnova[[i]][[length(tempAnova[[i]])]][[1]][numRow,"Mean Sq"] <- tempAnova[[i]][[length(tempAnova[[i]])]][[1]][numRow,"Sum Sq"]/dfError
                    tempAnova[[i]][[length(tempAnova[[i]])]][[1]][1:(numRow - 1),"F value"] <- tempAnova[[i]][[length(tempAnova[[i]])]][[1]][1:(numRow - 1),"Mean Sq"]/tempAnova[[i]][[length(tempAnova[[i]])]][[1]]["Residuals","Mean Sq"]
                    tempAnova[[i]][[length(tempAnova[[i]])]][[1]][1:(numRow - 1), "Pr(>F)"] <- pf(tempAnova[[i]][[length(tempAnova[[i]])]][[1]][1:(numRow - 1),"F value"], tempAnova[[i]][[length(tempAnova[[i]])]][[1]][1:(numRow - 1),"Df"], dfError, lower.tail = FALSE)
               } else {
                    tempAnova[[i]] <- summary(result)
                    tempAnova[[i]][[1]]["Df"] <- tempresult[[1]]["Df"]
                    tempAnova[[i]][[1]]["Mean Sq"] <- tempAnova[[i]][[1]]["Sum Sq"]/tempAnova[[i]][[1]]["Df"]
                    numEffects <- nrow(tempAnova[[i]][[1]])-1
                    dfError <- tempAnova[[i]][[1]][nrow(tempAnova[[i]][[1]]),"Df"]
                    tempAnova[[i]][[1]][1:numEffects, "F value"] <- tempAnova[[i]][[1]][1:numEffects,"Mean Sq"]/tempAnova[[i]][[1]]["Residuals","Mean Sq"]
                    tempAnova[[i]][[1]][1:numEffects, "Pr(>F)"] <- pf(tempAnova[[i]][[1]][1:numEffects,"F value"], tempAnova[[i]][[1]][1:numEffects,"Df"], dfError, lower.tail = FALSE)
               }
          } else { tempAnova[[i]] <- summary(result) }
          
          
          #availableDesign <- c("CRD", "RCBD", "LSD", "SplitCRD", "SplitRCBD", "SplitLSD", "Strip", "Split2CRD", "Split2RCBD", "Split2LSD", "Strip-Split",     "Split3CRD", "Split3RCBD", "Split3LSD", "Strip-Split2")
          if (design == "RCBD" || design == "SplitRCBD" || design == "Strip" || design == "Split2RCBD" || design == "Strip-Split" || design == "Split3RCBD" || design == "Strip-Split2") {
               if (attr(summary(result), "class")[[1]] == "summary.aovlist") { 
                    if (!is.na(match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]][[1]]))))) {
                         indexNum <- match(set, trimStrings(rownames(tempAnova[[i]][[1]][[1]])))
                         indexDen <- match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]][[1]])))
                         tempAnova[[i]][[1]][[1]][indexNum, 4] <- tempAnova[[i]][[1]][[1]][indexNum,3]/tempAnova[[i]][[1]][[1]][indexDen,3]
                         tempAnova[[i]][[1]][[1]][indexNum, 5] <- pf(tempAnova[[i]][[1]][[1]][indexNum,"F value"], tempAnova[[i]][[1]][[1]][indexNum,"Df"], tempAnova[[i]][[1]][[1]][indexDen,"Df"], lower.tail = FALSE)
                    }
               } else {
                    if (!is.na(match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]]))))) {
                         indexNum <- match(set, trimStrings(rownames(tempAnova[[i]][[1]])))
                         indexDen <- match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]])))
                         tempAnova[[i]][[1]][indexNum, 4] <- tempAnova[[i]][[1]][indexNum,3]/tempAnova[[i]][[1]][indexDen,3]
                         tempAnova[[i]][[1]][indexNum, 5] <- pf(tempAnova[[i]][[1]][indexNum,"F value"], tempAnova[[i]][[1]][indexNum,"Df"], tempAnova[[i]][[1]][indexDen,"Df"], lower.tail = FALSE)
                    }
               }
               
          }
          
          if (design == "LSD" || design == "SplitLSD" || design == "Split2LSD" || design == "Split3LSD") {
               if (attr(summary(result), "class")[[1]] == "summary.aovlist") { 
                    if (!is.na(match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]][[1]]))))) {
                         indexNum <- match(set, trimStrings(rownames(tempAnova[[i]][[1]][[1]])))
                         indexDen1 <- match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]][[1]])))
                         indexDen2 <- match(paste(set, rep2, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]][[1]])))
                         if (tempAnova[[i]][[1]][[1]][indexDen1,3] > tempAnova[[i]][[1]][[1]][indexDen2,3]) { indexDen <- indexDen1 
                         } else { indexDen <- indexDen2 }
                         tempAnova[[i]][[1]][[1]][indexNum, 4] <- tempAnova[[i]][[1]][[1]][indexNum,3]/tempAnova[[i]][[1]][[1]][indexDen,3]
                         tempAnova[[i]][[1]][[1]][indexNum, 5] <- pf(tempAnova[[i]][[1]][[1]][indexNum,"F value"], tempAnova[[i]][[1]][[1]][indexNum,"Df"], tempAnova[[i]][[1]][[1]][indexDen,"Df"], lower.tail = FALSE)
                    }
               } else {
                    if (!is.na(match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]]))))) {
                         indexNum <- match(set, trimStrings(rownames(tempAnova[[i]][[1]])))
                         indexDen1 <- match(paste(set, rep1, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]])))
                         indexDen2 <- match(paste(set, rep2, sep = ":"), trimStrings(rownames(tempAnova[[i]][[1]])))
                         if (tempAnova[[i]][[1]][indexDen1,3] > tempAnova[[i]][[1]][indexDen2,3]) { indexDen <- indexDen1 
                         } else { indexDen <- indexDen2 }
                         tempAnova[[i]][[1]][indexNum, 4] <- tempAnova[[i]][[1]][indexNum,3]/tempAnova[[i]][[1]][indexDen,3]
                         tempAnova[[i]][[1]][indexNum, 5] <- pf(tempAnova[[i]][[1]][indexNum,"F value"], tempAnova[[i]][[1]][indexNum,"Df"], tempAnova[[i]][[1]][indexDen,"Df"], lower.tail = FALSE)
                    }
               }
          }
          tempAOVTable <- ConstructAOVTable(tempAnova[[i]])
          
          # rename the anova table
          rownames(tempAOVTable) <- gsub("Error", "Pooled Error",trimStrings(rownames(tempAOVTable)))
          if (!is.null(rep1)) {
               index <- match(paste(c(set, rep1), collapse = ":", sep = ""), trimStrings(rownames(tempAOVTable)))
               if (!is.na(index)) {
                    rownames(tempAOVTable)[index] <- paste(rep1, "within", set)
                    tempAOVTable <- rbind(tempAOVTable[c(1,index),], tempAOVTable[-I(match(c(1, index), 1:nrow(tempAOVTable))),])
               }
               rm(index)
          }

          if (!is.null(rep2)) {
               index <- match(paste(c(set, rep2), collapse = ":", sep = ""), trimStrings(rownames(tempAOVTable)))
               if (!is.na(index)) {
                    rownames(tempAOVTable)[index] <- paste(rep2, "within", set)
                    tempAOVTable <- rbind(tempAOVTable[c(1,index),], tempAOVTable[-I(match(c(1, index), 1:nrow(tempAOVTable))),])
               }
               rm(index)
          }
          
          # -- CREATE THE RESIDUAL DATA AND PREDICTED VALUES -- #
          residNfittedData <- NULL
          residNfittedData <- data.frame(PredictedValues(result))
          if (inherits(result, what = "aovlist")) { residNfittedData <- data.frame(residNfittedData,proj(result)[[length(result)]][,"Residuals"])
          } else { residNfittedData <- data.frame(residNfittedData, residuals(result)) }
          colnames(residNfittedData) <- c(paste(respvar[i],"pred", sep = "_"), paste(respvar[i],"resid", sep = "_"))
          
          # -- CREATE THE DIAGNOSTIC PLOT -- #
          if (!is.null(outputPath)) {
               png(filename = paste(outputPath, design,"DiagPlot_", respvar[i], ".png", sep = ""), width = 480, height = 240)
               params <- par(mfrow = c(1,2), bg = "white")
               plot(residNfittedData[,(ncol(residNfittedData)-1)], residNfittedData[,ncol(residNfittedData)], main = paste("Residual vs Fitted:\n", respvar[i], sep = ""), xlab = "Fitted Values", ylab = "Residuals") 
               qqnorm(residNfittedData[,ncol(residNfittedData)])
               qqline(residNfittedData[,ncol(residNfittedData)])
               #plot(result)
               par(params)
               dev.off()
          }

          # -- PERFORM NORMALITY TEST AND/HOMOGENEITY OF VARIANCES -- #         
          if (normality || homogeneity) {
               assumpData <- data.frame(CombineFactorLevels(data = tempData, concatVar = allFactor, targetName = "factor")["factor"], residNfittedData[ncol(residNfittedData)])
               
               # --- PRINTING RESULTS OF TEST FOR HOMOGENEITY OF VARIANCES --- #
               if (homogeneity) {
                    capture.output(bartlett.result <- HeteroskedasticityTest(data = assumpData, var = paste(names(assumpData)[2]), grp = "factor", method = c("bartlett")))
                    cat("Bartlett's Test for Homogeneity of Variances\n")
                    printDataFrame(bartlett.result[,3:ncol(bartlett.result)])
                    cat("\n")
                    rm(bartlett.result)
               }
               
               # --- PRINTING RESULT OF SHAPIRO WILK TEST --- #
               if (normality) {
                    if (nrow(assumpData) >= 3 && nrow(assumpData) <= 5000) {
                         NormalityTest(data = assumpData, var = paste(names(assumpData)[2]), grp = NULL, method = c("swilk"))
                         cat("\n")
                    }
               }
               rm(assumpData)
          }
          
     
          # --- PRINTING OF ANOVA TABLE --- #
          if (is.null(contrastOption)) {
               cat("ANOVA TABLE\nResponse Variable: ", respvar[i], "\n", sep = "")
               printAOVTable(tempAOVTable)
               if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
               cat("\n")
          } else { 
               ContrastCompute(data = tempData, aovTable = tempAnova[[i]], mymodel, mymodel2,contrast.option = contrastOption)
               if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
               cat("\n")
          }
          
          # --- PRINTING OF SUMMARY STATISTICS --- #
          summaryStat <- NULL
          if (designChoice <= 3) {
               #if (designChoice == 1 && is.list(replications(formula, tempData))) { 
                    summaryTable <- suppressWarnings(model.tables(result, "means", se = FALSE))
               #} else { summaryTable <- suppressWarnings(model.tables(result, "means", se = TRUE)) }
               grandMean <- summaryTable$tables[[1]]
               summaryStat <- rbind(summaryStat, data.frame(((sqrt(tempAnova[[i]][[1]][nrow(tempAnova[[i]][[1]]),3])/grandMean) * 100))) 
               rownames(summaryStat)[nrow(summaryStat)] <- paste("CV(%)", sep = "")
               summaryStat <- t(rbind(summaryStat, grandMean))          	
          } else {
               grandMean <- mean(tempData[, respvar[i]], na.rm = TRUE)
               for (j in (1:length(tempAnova[[i]]))) { 
                    summaryStat <- rbind(summaryStat, data.frame(((sqrt(tempAnova[[i]][[j]][[1]][nrow(tempAnova[[i]][[j]][[1]]),3])/grandMean) * 100))); 
                    rownames(summaryStat)[nrow(summaryStat)] <- paste("CV(",letters[j],")%", sep = "")
               }
               summaryStat <- t(rbind(summaryStat, grandMean))
          }
          colnames(summaryStat)[ncol(summaryStat)] <- paste(respvar[i], "Mean")
          cat("Summary Statistics\n")
          printDataFrame(as.data.frame(summaryStat))
          if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
          cat("\n")
          
          # if (!estimatedData) {
          #      if (designChoice == 1 || designChoice == 2 || designChoice == 3) {
          #           if (!is.null(summaryTable$se)) {
          #                stdErrTable <- data.frame(Effects = names(unlist(summaryTable$se)),StdErr = unlist(summaryTable$se))
          #                rownames(stdErrTable) <- 1:nrow(stdErrTable)
          #                cat("Standard Errors\n")
          #                printDataFrame(stdErrTable)
          #                cat("\n")
          #           }
          #      }
          # }
          
          # --- DETERMINE THE EFFECTS WHICH ARE SIGNIFICANT --- #
          sigEffect <- SignificantEffect(tempAOVTable, alpha = sig)
          if (!is.null(sigEffect)) { 
               sigEffect <- trimStrings(sigEffect)
               deletedEffect <- FALSE
               if (!is.na(match(set, sigEffect))) {
                    sigEffect <- sigEffect[-I(match(set, sigEffect))]
                    deletedEffect <- TRUE
                    if (length(sigEffect) == 0) sigEffect <- NULL
               }
               if (!is.null(rep1) && !is.null(sigEffect)) {
                    if (!is.na(match(paste(rep1, "within", set), sigEffect))) {
                         sigEffect <- sigEffect[-I(match(paste(rep1, "within", set), sigEffect))]
                         deletedEffect <- TRUE
                         if (length(sigEffect) == 0) sigEffect <- NULL
                    }
               }
               if (!is.null(sigEffect) && !is.null(rep2)) {
                    if (!is.na(match(rep2, sigEffect))) {
                         sigEffect <- sigEffect[-I(match(rep2, sigEffect))]
                         deletedEffect <- TRUE
                         if (length(sigEffect) == 0) sigEffect <- NULL
                    }
               }
               if (!deletedEffect) rvWithSigEffect <- c(rvWithSigEffect, respvar[i])
          }
          
          # --- PRINT THE TABLE OF MEANS --- #
          if (is.null(sigEffect)) {
               cat("Table of Means\n")
               if (length(allFactor) == 1) {
                    tableMeans <- as.data.frame.table(summaryTable$tables[[length(summaryTable$tables)]])
                    colnames(tableMeans)[ncol(tableMeans)] <- paste(respvar[i]," Means", sep = "")
                    printDataFrame(tableMeans)
                    if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
               } else {
                    if (designChoice <= 3) {
                         print(ftable(summaryTable$tables[[length(summaryTable$tables)]]))
                         if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    } else {
                         factorInOrder <- unlist(lapply(tempData[allFactor], nlevels))[order(unlist(lapply(tempData[allFactor], nlevels)))]
                         tableMeans <-eval(parse(text = paste("ftable(tapply(tempData[,'",respvar[i],"'], list(tempData[,'", paste(names(factorInOrder), collapse = "'],tempData[,'", sep = ""),"']), mean))", sep = "")))
                         names(attr(tableMeans, "row.vars")) <- names(factorInOrder[1:(length(allFactor) - 1)])
                         names(attr(tableMeans, "col.vars")) <- names(factorInOrder[length(allFactor)])
                         print(tableMeans)
                         if (estimatedData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    }
               }
               cat("\n\n")
          } else {
               if (length(allFactor) > 1) {
                    highestInteraction <- paste(allFactor, collapse = ":", sep = "")
                    if (is.na(match(highestInteraction, sigEffect))) {
                         cat("Table of Means\n")
                         if (designChoice <= 3) { print(ftable(summaryTable$tables[[length(summaryTable$tables)]]))
                         } else {
                              factorInOrder <- unlist(lapply(tempData[allFactor], nlevels))[order(unlist(lapply(tempData[allFactor], nlevels)))]
                              tableMeans <-eval(parse(text = paste("ftable(tapply(tempData[,'",respvar[i],"'], list(tempData[,'", paste(names(factorInOrder), collapse = "'],tempData[,'", sep = ""),"']), mean))", sep = "")))
                              names(attr(tableMeans, "row.vars")) <- names(factorInOrder[1:(length(allFactor) - 1)])
                              names(attr(tableMeans, "col.vars")) <- names(factorInOrder[length(allFactor)])
                              print(tableMeans)
                         }
                         cat("\n\n")
                    }
               }
          } ## END IF ELSE STMT
          
          # --- PRINT PAIRWISE MEANCOMPARISON RESULT --- #
          if (!is.null(sigEffect)) {
               if (!is.na(match(respvar[i], pwVar))) {
                    for (j in (1:length(sigEffect))) {
                         pairwiseComparison(tempAnova[[i]], design, trimStrings(sigEffect[j]), data = tempData, respvar[i], pwTest, siglevel = sig)
                    }     
               } else {
                    for (j in (1:length(sigEffect))) {
                         pairwiseComparison(tempAnova[[i]], design, trimStrings(sigEffect[j]), data = tempData, respvar[i], pwTest = NULL, siglevel = sig)
                    }     
               }
          }
          pwOption[[i]] <- list(rv = respvar[i], test = pwTest, sigEffect = sigEffect)
          cat("\n")	
          
          tempNewData <- cbind(tempData, residNfittedData)
          # -- save the dataset
          if (estimatedData) {
               tempNewData <- merge(tempData[,-I(match(respvar[i],names(tempData)))], tempNewData)
               commonCol <- match(names(tempData),names(tempNewData))
               tempData <- cbind(tempNewData[,commonCol], tempNewData[,-I(commonCol)])
          } else { if (!is.null(tempNewData)) tempData <- merge(tempData, tempNewData) }
          
     } ### end stmt --- for (i in (1:length(respvar)))
     
     options(show.signif.stars = prev.option)
     return(invisible(list(data = tempData, aovObject = aovresult, rvWithSigEffect = rvWithSigEffect, aovTable = tempAnova, pwOption = pwOption, model = modelRHS, model2 = modelRHS2, alpha = sig)))
     
     
}