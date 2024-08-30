# ----------------------------------------------------------------------
# pairwiseAmong: Function for displaying the mean comparison
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 
# ----------------------------------------------------------------------

pairwiseAmong <- function(data, respvar, typeTest, trmt, dfError, MSError, siglevel) UseMethod("pairwiseAmong")

pairwiseAmong.default <- function(data, respvar, typeTest, trmt, dfError, MSError, siglevel) {
     options(width = 5000)
     if (nlevels(data[,trmt]) > 26) {
          command <- paste(typeTest, ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = FALSE, pwOrder = 'trmt')", sep = "")
          eval(parse(text = command))
     } else {
          if (length(typeTest) == 1) {
               command <- paste(typeTest, ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
               eval(parse(text = command))
          } else {
               singleVal <- na.omit(match(typeTest, c("LSD", "HSD", "scheffe")))
               if (length(singleVal) != 0) {
                    summaryTable <- NULL
                    summaryStat <- NULL
                    for (i in (1:length(singleVal))) {
                         command <- paste(c("LSD", "HSD", "scheffe")[singleVal[i]], ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
                         if (length(singleVal) == 1) { eval(parse(text = command))
                         } else {
                              capture.output(resultPW <- eval(parse(text = command)))
                              if (is.null(summaryTable)) {
                                   summaryTable <- resultPW$summary[,c(1:3,5)]
                                   colnames(summaryTable)[4] <- c("LSD", "HSD", "scheffe")[singleVal[i]]
                                   summaryStat <- data.frame(rbind(resultPW$tabValue, resultPW$testStat))
                                   colnames(summaryStat) <- c("LSD", "HSD", "scheffe")[singleVal[i]]
                                   rownames(summaryStat) <- c("Critical Value", "Test Statistic")
                              } else {
                                   summaryTable <- cbind(summaryTable, resultPW$summary[,5])
                                   colnames(summaryTable)[3+i] <- c("LSD", "HSD", "scheffe")[singleVal[i]]
                                   summaryStat <- cbind(summaryStat, data.frame(rbind(resultPW$tabValue, resultPW$testStat)))
                                   colnames(summaryStat)[i] <- c("LSD", "HSD", "scheffe")[singleVal[i]]
                              }
                         }
                    }
                    
                    if (length(singleVal) > 1) {
                         maxWidthEntry <- max(nchar(dfError), nchar(round(MSError,0))) + 7
                         cat("\n")
                         cat(formatC("Alpha", format = "s", width = 25, flag = "-"), formatC(siglevel, format = "f", digits = 2, width = maxWidthEntry , flag = "#"), "\n", sep = "")
                         cat(formatC("Error Degrees of Freedom", format = "s", width = 25, flag = "-"), formatC(dfError, format = "d", width = maxWidthEntry , flag = "#"), "\n", sep = "")
                         cat(formatC("Error Mean Square", format = "s", width = 25, flag = "-"), formatC(MSError, format = "f", digits = 4, width = maxWidthEntry , flag = "#"), "\n\n", sep = "")
                         
                         maxWidthEntry <- nchar(max(round(summaryStat,0))) + 7
                         maxWidthLabel <- max(nchar(rownames(summaryStat))) + 2
                         
                         cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                         cat(formatC("", format = "s", width = maxWidthLabel, flag = "-"), formatC(colnames(summaryStat), format = "s", width = maxWidthEntry, flag = "#"), "\n", sep = "")
                         cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                         for (i in (1:2)) {
                              cat(formatC(rownames(summaryStat)[i], format = "s", width = maxWidthLabel, flag = "-"), sep = "")
                              for (j in (1:ncol(summaryStat))) {
                                   cat(formatC(summaryStat[i,j], format = "f", digits = 4, width = maxWidthEntry, flag = "#"),sep = "")
                              }
                              cat("\n")
                         }
                         cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                         
                         cat("\nSummary: \n")
                         printDataFrame(summaryTable, digits = 4)
                         cat("Means with the same letter are not significantly different\n\n")
                    } 
                    
               } 
               
               rangeVal <- na.omit(match(typeTest, c("duncan", "SNK")))
               if (length(rangeVal) != 0) {
                    for (i in (1:length(rangeVal))) {
                         command <- paste(c("duncan", "SNK")[rangeVal[i]], ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
                         eval(parse(text = command))
                    }
               }
          }
     }
     
} ## END FUNCTION
