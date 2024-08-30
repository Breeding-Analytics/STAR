# --------------------------------------------------------------------------------------------------
# printBIBPairwise
# Print the result of the pairwise mean comparison for Balanced Incomplete Block Design
# Date Created: 10.17.2013
# Created by: AAG for International Rice Research Institute
# --------------------------------------------------------------------------------------------------

printBIBDPairwise <- function(object, dfError, MSError, method, trmtLabel = NULL) UseMethod("printBIBDPairwise")
     
printBIBDPairwise.default <- function(object, dfError, MSError, method, trmtLabel = NULL) {
     cat("Pairwise Mean Comparison")
     if (!is.null(trmtLabel)) {  cat(" of",trmtLabel) }
     cat("\n\n")
     
     for (i in (1:length(method))) {
          index <- match(method[i], names(object))
          cat(object[[index]]$method, "\n\n")
          cat(formatC("Alpha", format = "s", width = 24, flag = "-"), formatC(object[[index]]$alpha, format = "f", digits = 2, width = 10, flag = "#"), "\n", sep = "")
          cat(formatC("Error Degrees of Freedom", format = "s", width = 24, flag = "-"), formatC(dfError, format = "d", width = 10, flag = "#"), "\n", sep = "")
          cat(formatC("Error Mean Square", format = "s", width = 24, flag = "-"), formatC(MSError, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
          if (names(object)[index] == "LSD" || names(object)[index] == "HSD") {
               cat(formatC("Critical Value", format = "s", width = 24, flag = "-"), formatC(object[[index]]$cval/object[[index]]$stderr, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
               cat(formatC(paste("Test Statistics"), format = "s", width = 24, flag = "-"), formatC(object[[index]]$cval, format = "f", digits = 4, width = 10, flag = "#"), "\n\n", sep = "")
          } else {
               cat("\n")
               cat("\n")
               tempTable <- cbind(c("Critical Value", "Test Statistic"),as.data.frame(object[[index]]$crange))
               as.data.frame(object[[index]]$crange)
               names(tempTable)[1] <- "Number of Means"
               printDataFrame(tempTable, digits = 4)
               cat("\n")
               #cat(formatC("Number of Means", format = "s", width = 24, flag = "-"), formatC(colnames(object[[index]]$crange), format = "d", width = 10, flag = "#"), "\n") 
               #cat(formatC("Critical Value", format = "s", width = 24, flag = "-"), formatC(object[[index]]$crange[1,], format = "f", digits = 4, width = 10, flag = "#"), "\n")
               #cat(formatC(paste("Test Statistics"), format = "s", width = 24, flag = "-"), formatC(object[[index]]$crange[2,], format = "f", digits = 4, width = 10, flag = "#"), "\n\n")
          }
          printDataFrame(object[[index]]$pw, border = TRUE)
          cat("Means with the same letter are not significantly different. \n\n\n")
     }
}