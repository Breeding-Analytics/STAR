# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# aovTest: Functions for performing ANOVA
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.19.2013
# Note: Include remarks when data is unbalanced, when no balanced data can be generated
# -------------------------------------------------------------------------------

AOVTest <- function(data, design, respvar, factor1, factor2 = NULL, factor3 = NULL, factor4 = NULL, rep1 = NULL, rep2 = NULL, set = NULL, descriptive = FALSE, normality = FALSE, homogeneity = FALSE, pwTest = NULL, pwVar = NULL, contrastOption = NULL, sig = 0.05, outputPath = NULL) UseMethod("AOVTest")

AOVTest.default <- function(data, design, respvar, factor1, factor2 = NULL, factor3 = NULL, factor4 = NULL, rep1 = NULL, rep2 = NULL, set = NULL, descriptive = FALSE, normality = FALSE, homogeneity = FALSE, pwTest = NULL, pwVar = NULL, contrastOption = NULL, sig = 0.05, outputPath = NULL) {
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
                      "Strip Plot Design", "Split-Split Plot in Completely Randomized Design",
                      "Split-Split Plot in Randomized Complete Block Design", 
                      "Split-Split Plot in Latin Square Design",
                      "Strip-Split Plot Design", 
                      "Split-Split-Split Plot in Completely Randomized Design", 
                      "Split-Split-Split Plot in Randomized Complete Block Design", 
                      "Split-Split-Split Plot in Latin Square Design", "Strip-Split-Split Plot Design")
     
     switch(designChoice,
{numfactor <- 1; numblk <- 0; modelRHS <- paste(paste(factor1, collapse = "*", sep = "")); modelRHS2 <- modelRHS},
{numfactor <- 1; numblk <- 1; modelRHS <- paste(rep1, " + ", paste(factor1, collapse = "*", sep = ""), sep = ""); modelRHS2 <- modelRHS},
{numfactor <- 1; numblk <- 2; modelRHS <- paste(rep1, " + ", rep2, " + ", paste(factor1, collapse = "*", sep = ""), sep = ""); modelRHS2 <- modelRHS},
{numfactor <- 2; numblk <- 1; modelRHS <- paste(paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = "") , ")/",  rep1, ")", sep = "");
 modelRHS2 <- paste(paste(c(factor1, factor2), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), sep = "")},
{numfactor <- 2; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error(", paste(c(rep1, factor1), collapse = ":", sep = ""), "/(", paste(factor1, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), sep = "")},
{numfactor <- 2; numblk <- 2; modelRHS <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error(", paste(c(rep1,rep2, factor1), collapse = ":", sep = ""), "/(", paste(factor1, collapse = "*", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + ", paste(c(rep1,rep2, factor1), collapse = ":", sep = ""), sep = "")},
{numfactor <- 2; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error((", paste(c(rep1,factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + ", paste(c(rep1,factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor2), collapse = ":", sep = ""), sep = "")},
{numfactor <- 3; numblk <- 1; modelRHS <- paste(paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), sep = "")},
{numfactor <- 3; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), sep = "")},
{numfactor <- 3; numblk <- 2; modelRHS <- paste(rep1 ," + ", rep2 ," + ",paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", rep2 ," + ",paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + ", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), sep = "")},
{numfactor <- 3; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(c(factor1, factor2), collapse = ":", sep = ""),"))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor2), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), sep = "")},
{numfactor <- 4; numblk <- 1; modelRHS <- paste(paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(factor3, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), sep = "")},
{numfactor <- 4; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""),         "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""),         "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(factor3, collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), sep = "")},
{numfactor <- 4; numblk <- 2; modelRHS <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(factor3, collapse = ":", sep = ""),"))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + ", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), " + ", paste(c(rep1, rep2, factor1, factor2, factor3), collapse = ":", sep = ""), sep = "")},
{numfactor <- 4; numblk <- 1; modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(c(factor1, factor2), collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(c(factor3), collapse = ":", sep = ""), "))", sep = "");
 modelRHS2 <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + ", paste(c(rep1, factor1), collapse = ":", sep = ""), " + ", paste(c(rep1, factor2), collapse = ":", sep = ""), " + ", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), sep = "")}
     )	
     
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
     #for (i in (1:length(allFactor))) { data[,allFactor[i]]	<- factor(data[,allFactor[i]])	 }
     for (i in (1:length(allFactor))) { 
          if (is.numeric(data[,allFactor[i]])) { data[,allFactor[i]] <- factor(data[,allFactor[i]])
          } else {
               if (suppressWarnings(all(!is.na(as.numeric(as.character(factor(trimStrings(data[,allFactor[i]])))))))) {
                    data[,allFactor[i]] <- factor(as.numeric(data[,allFactor[i]]))
               } else { data[,allFactor[i]] <- factor(trimStrings(data[,allFactor[i]])) }
          }
     }
     
     if (!is.null(rep1)) { data[,rep1] <- factor(data[,rep1]) }
     if (!is.null(rep2)) { data[,rep2] <- factor(data[,rep2]) }
     
     if (is.null(set)) {
          set <- make.unique(c(names(data), "setVar"))[length(make.unique(c(names(data), "setVar")))]
          data <- cbind(data, setVar = "1")
          addSet <- TRUE
     } else {
          addSet <- FALSE   
          data[,set] <- factor(data[,set])
     }
     tempAllData <- data
     
     cat("Analysis of Variance","\n",designTitle[designChoice],"\n\n", sep = "")
     
     for (i in (1:length(respvar))) {
          width1 <- 37 + nchar(respvar[i])
          cat(paste(rep("=", width1), collapse = ""), "\n")
          cat("ANALYSIS FOR RESPONSE VARIABLE:", respvar[i],"\n")
          cat(paste(rep("=", width1), collapse = ""), "\n\n")
          
          aovresult[[i]] <- list()
          tempAnova[[i]] <- list()          
          pwOption[[i]] <- list()
          estimateAllData <- FALSE
          tempNewData <- NULL
          
          for (z in (1:nlevels(data[,match(set, names(data))]))) {
               proceedPrinting <- TRUE
               estimateData <- FALSE
               setLabel <- levels(data[,match(set, names(data))])[z]
               if (!addSet) {
                    width2 <- 5 + nchar(set) + nchar(setLabel)
                    cat(paste(rep("-", width2), collapse = ""), "\n")
                    cat(set, "=", setLabel, "\n")
                    cat(paste(rep("-", width2), collapse = ""), "\n\n")
                    rm(width2)
               }
               
               origData <- tempAllData[tempAllData[,set] == levels(tempAllData[,match(set, names(tempAllData))])[z],]
               
               # determine if the data set contains an observation
               if (length(na.omit(origData[,respvar[i]])) == 0) {
                    cat("Error: The data set does not contain any observations.\n\n")
                    next
               }
               
               # converts to a factor
               for (y in (1:length(allFactor))) { origData[,allFactor[y]]     <- factor(origData[,allFactor[y]])	 }
               if (!is.null(rep1)) { origData[,rep1] <- factor(origData[,rep1]) }
               if (!is.null(rep2)) { origData[,rep2] <- factor(origData[,rep2]) }
               tempData <- origData
               
               formula <- paste(respvar[i]," ~ ", paste(c(rep1, rep2, allFactor), collapse = ":", sep = ""), sep = "")
               
               # determine if balanced data
               if (designChoice != 1) {
                    rawData <- tempData[,sort(match(c(respvar[i], allFactor, rep1, rep2), names(tempData)))]
                    if (nrow(na.omit(rawData)) == 0) {
                         cat("ERROR: Data set has no observation\n\n", sep = "")
                         next
                    }
                    if (designChoice != 3 && designChoice != 6 && designChoice != 10 && designChoice != 14) {
                         if (is.list(replications(formula, rawData))) {
                              if (max(replications(formula, rawData)[[1]]) != 1) {
                                   cat("ERROR: Cannot perform ANOVA for balanced data. The data set for response variable '", respvar[i],"' is unbalanced.\n\n", sep = "")
                                   next
                              }
                              newData <- GenerateBalanceData(tempData, respvar[i], allFactor, c(rep1, rep2) , design)
                              if ((1 - (length(na.omit(newData[,respvar[i]]))/nrow(newData))) > 0.10) { 
                                   cat("ERROR: Cannot perform ANOVA for balanced data for response variable '", respvar[i], "'. Too many missing values.\n\n", sep = "")
                                   next 
                              }
                              tempData <- estMissData(design, data = newData, respvar[i], factor1, factor2, factor3, factor4, rep1, rep2)
                              estimateData <- TRUE
                              estimateAllData <- TRUE
                         } else { estimateData <- FALSE }     
                    } else {
                         # for LSD design structure
                         formula1 <- paste(respvar[i]," ~ ", paste(c(rep1, allFactor), collapse = ":", sep = ""), sep = "") 
                         if (is.list(replications(formula1, rawData))) {
                              if (max(replications(formula1, rawData)[[1]]) != 1) {
                                   cat("ERROR: Cannot perform ANOVA for balanced data. The data set for response variable '", respvar[i],"' is unbalanced.\n\n", sep = "")
                                   next
                              }
                              newData <- GenerateBalanceData(tempData, respvar[i], allFactor, c(rep1, rep2) , design)
                              if ((1 - (length(na.omit(newData[,respvar[i]]))/nrow(newData))) > 0.10) { 
                                   cat("ERROR: Cannot perform ANOVA for balanced data for response variable '", respvar[i], "'. Too many missing values.\n\n", sep = "")
                                   next 
                              }
                              tempData <- estMissData(design, data = newData, respvar[i], factor1, factor2, factor3, factor4, rep1, rep2)
                              estimateData <- TRUE
                              estimateAllData <- TRUE
                         } else { estimateData <- FALSE }
                    }
               } else { 
                    rawData <- tempData[,sort(match(c(respvar[i], allFactor, rep1, rep2), names(tempData)))]
                    if (nrow(na.omit(rawData)) == 0) {
                         cat("ERROR: Data set has no observation\n\n", sep = "")
                         next
                    }
                    estimateData <- FALSE 
               }
               
               # determine if there is at least two observation per treatment levels or treatment combination 
               
               formulaNew <- paste(respvar[i]," ~ ", paste(allFactor, collapse = ":", sep = ""), sep = "")
               if (design != "CRD") {
                    if (replications(formulaNew, tempData[,sort(match(c(respvar[i], allFactor), names(tempData)))]) == 1) {
                         cat("Error: There is only one observation per treatment level.\n\n", sep = "")
                         next
                    } 
               }
               
               
               # determine if the variable is constant
               if (var(tempData[,respvar[i]]) == 0) {
                    cat("ERROR: The data for the response variable '", respvar[i], "' is constant.\n\n", sep = "")
                    next 
               }
               
               # define the model
               modelLHS <- paste(respvar[i], "~")
               mymodel <- paste(modelLHS, modelRHS)
               mymodel2 <- paste(modelLHS, modelRHS2)
               if (estimateData) { tempresult <- summary(suppressWarnings(aov(formula(mymodel), data = origData)))     }
               resultTry <- try(result <- suppressWarnings(aov(formula(mymodel), tempData)), silent = TRUE)
               if(all(class(resultTry) == "try-error")) {
                    msg <- trimStrings(strsplit(resultTry, ":")[[1]])
                    msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
                    cat("Error: ", msg, "\n\n", sep = "")
                    next
               }
               
               aovresult[[i]][[z]] <- result
               
               if (estimateData) {
                    if (attr(summary(result), "class")[[1]] == "summary.aovlist") {
                         tempAnova[[i]][[z]] <- summary(result) 
                         numRow <- nrow(tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]])
                         dfError <- tempresult[[length(tempresult)]][[1]][nrow(tempresult[[length(tempresult)]][[1]]),"Df"]
                         tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][numRow,"Df"] <- dfError
                         tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][numRow,"Mean Sq"] <- tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][numRow,"Sum Sq"]/dfError
                         tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][1:(numRow - 1),"F value"] <- tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][1:(numRow - 1),"Mean Sq"]/tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]]["Residuals","Mean Sq"]
                         tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][1:(numRow - 1), "Pr(>F)"] <- pf(tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][1:(numRow - 1),"F value"], tempAnova[[i]][[z]][[length(tempAnova[[i]][[z]])]][[1]][1:(numRow - 1),"Df"], dfError, lower.tail = FALSE)
                    } else {
                         tempAnova[[i]][[z]] <- summary(result)
                         tempAnova[[i]][[z]][[1]]["Df"] <- tempresult[[1]]["Df"]
                         tempAnova[[i]][[z]][[1]]["Mean Sq"] <- tempAnova[[i]][[z]][[1]]["Sum Sq"]/tempAnova[[i]][[z]][[1]]["Df"]
                         numEffects <- nrow(tempAnova[[i]][[z]][[1]])-1
                         dfError <- tempAnova[[i]][[z]][[1]][nrow(tempAnova[[i]][[z]][[1]]),"Df"]
                         tempAnova[[i]][[z]][[1]][1:numEffects, "F value"] <- tempAnova[[i]][[z]][[1]][1:numEffects,"Mean Sq"]/tempAnova[[i]][[z]][[1]]["Residuals","Mean Sq"]
                         tempAnova[[i]][[z]][[1]][1:numEffects, "Pr(>F)"] <- pf(tempAnova[[i]][[z]][[1]][1:numEffects,"F value"], tempAnova[[i]][[z]][[1]][1:numEffects,"Df"], dfError, lower.tail = FALSE)
                    }
               } else { tempAnova[[i]][[z]] <- summary(result) }
               tempAOVTable <- ConstructAOVTable(tempAnova[[i]][[z]])
               
               # -- PRINTING CLASS LEVEL INFORMATION -- #
               #ClassInformation(tempData[, c(factor1, factor2, factor3, factor4, rep1, rep2, respvar[i])], respvar = respvar[i])
               ClassInformation(origData[, c(factor1, factor2, factor3, factor4, rep1, rep2, respvar[i])], respvar = respvar[i])
               cat("\n\n")
               
               # --- PRINTING DESCRIPTIVE STATISTICS --- #     
               if (descriptive) { 
                    DescriptiveStatistics(data = tempData, var = respvar[i], grp = NULL, statistics = c("n", "mean", "sd", "min", "max"))
                    if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    cat("\n")
               }
               
               # -- CREATE THE RESIDUAL DATA AND PREDICTED VALUES -- #
               residNfittedData <- NULL
               residNfittedData <- data.frame(PredictedValues(result))
               if (inherits(result, what = "aovlist")) { residNfittedData <- data.frame(residNfittedData,proj(result)[[length(result)]][,"Residuals"])
               } else { residNfittedData <- data.frame(residNfittedData, residuals(result)) }
               colnames(residNfittedData) <- c(paste(respvar[i],"pred", sep = "_"), paste(respvar[i],"resid", sep = "_"))
               
               # -- CREATE THE DIAGNOSTIC PLOT -- #
               if (!is.null(outputPath)) {
                    if (addSet) { png(filename = paste(outputPath, design,"DiagPlot_", respvar[i], ".png", sep = ""), width = 480, height = 240)
                    } else {
                         tempSetLabel <- setLabel
                         splitSetLabel <- strsplit(setLabel, split = "")
                         if(!is.na(match("/", splitSetLabel[[1]]))) {
                              splitSetLabel <- gsub("/", "_", splitSetLabel[[1]])
                              tempSetLabel <- capture.output(cat(splitSetLabel, sep = ""))
                         }
                         png(filename = paste(outputPath, design,"DiagPlot_", tempSetLabel,"_", respvar[i], ".png", sep = ""), width = 480, height = 240) 
                    }
                    params <- par(mfrow = c(1,2), bg = "white")
                    if (addSet) { plot(residNfittedData[,(ncol(residNfittedData)-1)], residNfittedData[,ncol(residNfittedData)], main = paste("Residual vs Fitted:\n", respvar[i], sep = ""), xlab = "Fitted Values", ylab = "Residuals") 
                    } else { plot(residNfittedData[,(ncol(residNfittedData)-1)], residNfittedData[,ncol(residNfittedData)], main = paste("Residual vs Fitted:\n", respvar[i], " for ", setLabel, sep = ""), xlab = "Fitted Values", ylab = "Residuals") }
                    qqnorm(residNfittedData[,ncol(residNfittedData)])
                    qqline(residNfittedData[,ncol(residNfittedData)])
                    par(params)
                    # params <- par(mfrow = c(2,2), bg = "white")
                    # plot(result)
                    # par(params)
                    dev.off()
               }
               
               # -- PERFORM NORMALITY TEST AND/HOMOGENEITY OF VARIANCES -- #
               if (homogeneity || normality) {
                    tempData2 <- tempData[!is.na(tempData[,respvar[i]]),]
                    assumpData <- data.frame(CombineFactorLevels(data = tempData2, concatVar = allFactor, targetName = "factor")["factor"], residNfittedData[ncol(residNfittedData)])
                    
                    # --- PRINTING RESULTS OF TEST FOR HOMOGENEITY OF VARIANCES --- #
                    if (homogeneity) {
                         capture.output(bartlett.result <- HeteroskedasticityTest(data = assumpData, var = paste(names(assumpData)[2]), grp = "factor", method = c("bartlett")))
                         cat("Test for Homogeneity of Variances\n")
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
               
               # determine if the column Mean Square contains a value less than 1.0e-15
               # then determine if the error term is 1.0e-15
               if (any(na.omit(tempAOVTable[,"Mean Sq"] <= 1.0e-15))) {
                    # determine if the ANOVA contains 2 or more error terms
                    if (attr(summary(result), "class")[[1]] == "summary.aovlist") {
                         for (j in (1:length(summary(result)))) {
                              if (j == 1) { startIndex <- 1 } else { startIndex <- errorIndex + 1 }
                              errorIndex <- match(paste("Error(",letters[j],")", sep = ""), rownames(tempAOVTable))
                              if (tempAOVTable[errorIndex,"Mean Sq"] <= 1.0e-15) {
                                   tempAOVTable[errorIndex, "Mean Sq"] <- 0   
                                   tempAOVTable[startIndex:(errorIndex-1), "Pr(>F)"] <- NaN
                                   tempAOVTable[startIndex:(errorIndex-1), "F value"] <- NaN
                                   proceedPrinting <- FALSE
                              }
                         }
                    } else {
                         if (tempAOVTable[match("Error", rownames(tempAOVTable)),"Mean Sq"] <= 1.0e-15) {
                              errorIndex <- match("Error", rownames(tempAOVTable))
                              tempAOVTable[errorIndex, "Mean Sq"] <- 0   
                              numEffects <- nrow(tempAOVTable) - 2
                              tempAOVTable[1:numEffects, "Pr(>F)"] <- NaN
                              tempAOVTable[1:numEffects, "F value"] <- NaN
                              proceedPrinting <- FALSE
                         }
                    }
               }
               
               # --- PRINTING OF ANOVA TABLE --- #
               if (is.null(contrastOption)) {
                    cat("ANOVA TABLE\nResponse Variable: ", respvar[i], "\n", sep = "")
                    printAOVTable(tempAOVTable)
                    if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    cat("\n")
               } else { 
                    ContrastCompute(data = tempData, aovTable = tempAnova[[i]][[z]], mymodel, mymodel2,contrast.option = contrastOption)
                    if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    cat("\n")
               }
               
               if (!proceedPrinting) { next }
               
               
               # --- PRINTING OF SUMMARY STATISTICS --- #
               summaryStat <- NULL
               if (designChoice <= 3) {
                    if (designChoice == 1 && is.list(replications(formula, tempData))) { 
                         summaryTable <- suppressWarnings(model.tables(result, "means", se = FALSE))
                    } else { summaryTable <- suppressWarnings(model.tables(result, "means", se = TRUE)) }
                    grandMean <- summaryTable$tables[[1]]
                    summaryStat <- rbind(summaryStat, data.frame(((sqrt(tempAnova[[i]][[z]][[1]][nrow(tempAnova[[i]][[z]][[1]]),3])/grandMean) * 100))) 
                    rownames(summaryStat)[nrow(summaryStat)] <- paste("CV(%)", sep = "")
                    summaryStat <- t(rbind(summaryStat, grandMean))               
               } else {
                    grandMean <- mean(tempData[, respvar[i]], na.rm = TRUE)
                    for (j in (1:length(tempAnova[[i]][[z]]))) { 
                         summaryStat <- rbind(summaryStat, data.frame(((sqrt(tempAnova[[i]][[z]][[j]][[1]][nrow(tempAnova[[i]][[z]][[j]][[1]]),3])/grandMean) * 100))); 
                         rownames(summaryStat)[nrow(summaryStat)] <- paste("CV(",letters[j],")%", sep = "")
                    }
                    summaryStat <- t(rbind(summaryStat, grandMean))
               }
               colnames(summaryStat)[ncol(summaryStat)] <- paste(respvar[i], "Mean")
               cat("Summary Statistics\n")
               printDataFrame(as.data.frame(summaryStat))
               if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
               cat("\n")
               
               if (!estimateData) {
                    if (designChoice == 1 || designChoice == 2 || designChoice == 3) {
                         if (!is.null(summaryTable$se)) {
                              stdErrTable <- data.frame(Effects = names(unlist(summaryTable$se)),StdErr = unlist(summaryTable$se))
                              rownames(stdErrTable) <- 1:nrow(stdErrTable)
                              cat("Standard Errors\n")
                              printDataFrame(stdErrTable)
                              cat("\n")
                         }
                    }
               }
               
               # --- DETERMINE THE EFFECTS WHICH ARE SIGNIFICANT --- #
               sigEffect <- SignificantEffect(tempAOVTable, alpha = sig)
               if (!is.null(sigEffect)) { 
                    sigEffect <- trimStrings(sigEffect)
                    deletedEffect <- FALSE
                    if (!is.null(rep1)) {
                         if (!is.na(match(rep1, sigEffect))) {
                              sigEffect <- sigEffect[-I(match(rep1, sigEffect))]
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
                    if (!deletedEffect) rvWithSigEffect <- unique(c(rvWithSigEffect, respvar[i]))
               }
               
               # --- PRINT THE TABLE OF MEANS --- #
               if (is.null(sigEffect)) {
                    cat("Table of Means\n")
                    if (length(allFactor) == 1) {
                         tableMeans <- as.data.frame.table(summaryTable$tables[[length(summaryTable$tables)]])
                         colnames(tableMeans)[ncol(tableMeans)] <- paste(respvar[i]," Means", sep = "")
                         printDataFrame(tableMeans)
                         if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                    } else {
                         if (designChoice <= 3) {
                              print(ftable(summaryTable$tables[[length(summaryTable$tables)]]))
                              if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
                         } else {
                              factorInOrder <- unlist(lapply(tempData[allFactor], nlevels))[order(unlist(lapply(tempData[allFactor], nlevels)))]
                              tableMeans <-eval(parse(text = paste("ftable(tapply(tempData[,'",respvar[i],"'], list(tempData[,'", paste(names(factorInOrder), collapse = "'],tempData[,'", sep = ""),"']), mean))", sep = "")))
                              names(attr(tableMeans, "row.vars")) <- names(factorInOrder[1:(length(allFactor) - 1)])
                              names(attr(tableMeans, "col.vars")) <- names(factorInOrder[length(allFactor)])
                              print(tableMeans)
                              if (estimateData) { cat("REMARK: Raw data and estimates of the missing values are used.\n") } 
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
                              pairwiseComparison(tempAnova[[i]][[z]], design, trimStrings(sigEffect[j]), data = tempData, respvar[i], pwTest, siglevel = sig)
                         }     
                    } else {
                         for (j in (1:length(sigEffect))) {
                              pairwiseComparison(tempAnova[[i]][[z]], design, trimStrings(sigEffect[j]), data = tempData, respvar[i], pwTest = NULL, siglevel = sig)
                         }     
                    }
               }
               
               #if (!is.null(sigEffect)) { 
               if (addSet) { pwOption[[i]][[z]] <- list(rv = respvar[i], test = pwTest, sigEffect = sigEffect, setvar = "NULL", set = levels(data[,match(set, names(data))])[z]) 
               } else { pwOption[[i]][[z]] <- list(rv = respvar[i], test = pwTest, sigEffect = sigEffect, setvar = set, set = levels(data[,match(set, names(data))])[z])  }
               #} else { pwOption[[i]][[z]] <- list(rv = respvar[i], test = pwTest, sigEffect = sigEffect) }
               cat("\n")
               
               # -- COMBINE DATASET -- #
               #if (z == 1) { tempNewData <- cbind(tempData, residNfittedData) 
               #} else { tempNewData <- rbind(tempNewData, cbind(tempData, residNfittedData)) }
               if (design == "CRD") {
                    tempData1 <- matrix(NA, nrow = nrow(tempData[is.na(tempData[,respvar[i]]),]), ncol = 2)
                    colnames(tempData1) <- colnames(residNfittedData)
                    tempNewData <- rbind(tempNewData, 
                                         rbind(cbind(tempData[!is.na(tempData[,respvar[i]]),], residNfittedData),
                                               cbind(tempData[is.na(tempData[,respvar[i]]),], tempData1)))
               } else { 
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
               }
               
               
          } ## end stmt -- for (z in (1:nlevels(data[,match(set, names(data))])))    
          
          # -- save the dataset
          if (estimateAllData) {
               tempNewData <- merge(tempAllData[,-I(match(respvar[i],names(tempAllData)))], tempNewData, all = TRUE)
               commonCol <- match(names(tempAllData),names(tempNewData))
               tempAllData <- cbind(tempNewData[,commonCol], tempNewData[,-I(commonCol)])
          } else { if (!is.null(tempNewData)) tempAllData <- merge(tempAllData, tempNewData, all = TRUE) }
          
          
     } ### end stmt --- for (i in (1:length(respvar)))
     
     options(show.signif.stars = prev.option)
     return(invisible(list(data = tempAllData, aovObject = aovresult, rvWithSigEffect = rvWithSigEffect, aovTable = tempAnova, pwOption = pwOption, model = modelRHS, model2 = modelRHS2, alpha = sig)))
     
} ### end stmt -- aovTest
