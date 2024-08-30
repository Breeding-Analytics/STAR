#----------------------------------------------------#
#Discriminant - Function for performing Discriminant Analysis
#
# data - name of R dataframe
# outputPath - folder where output(s) will be saved
# vars - vector of names of the numeric variables
# grpVar - name of the grouping variable
# descriptive - logical; whether descriptive statistics will be displayed or not
# doNormalTest - logical; whether multivariate normality test will be performed or not
# doBoxM - logical; whether Box's M test will be performed or not
# priorProb - vector of prior probabilities to be used, if not taken from the observed freq from the data
# classifyNew - logical; whether  observations in a new file will be classified or not
# data2 - name of source file of observations to be classified
#----------------------------------------------------#

Discriminant <- function(data, outputPath, vars, grpVar, descriptive = FALSE, 
                        doNormalTest = FALSE, doBoxM = FALSE, priorProb = NULL,  
                        classifyNew = FALSE, data2)
   UseMethod("Discriminant")

Discriminant.default <- function(data, outputPath, vars, grpVar, descriptive = FALSE, 
                        doNormalTest = FALSE, doBoxM = FALSE, priorProb = NULL,  
                        classifyNew = FALSE, data2) {
  
  #reads data
  if (is.character(data)) {
    nameData <- data
    data <- eval(parse(text = data))
  } else { 
    nameData <- paste(deparse(substitute(data))) 
    tempData <- data
  }
    
  #converts to factor the grouping variable
  tempData[,grpVar] <- factor(tempData[,grpVar]) 
  
  # -- PRINTING CLASS LEVEL INFORMATION -- #
  ClassInformation(tempData[, c(grpVar, vars)], respvar = c(grpVar, vars))
  cat("\n")
  
  data <- na.omit(data[,c(vars,grpVar)])
  if (nrow(data) == 0) { stop("Data does not contain any observations.") }
  
  #converts to factor the grouping variable
  data[,grpVar] <- factor(data[,grpVar]) 
 
  results <- list()
  
  ###OPTIONS
  #DESCRIPTIVE     
  if (descriptive) { 
    DescriptiveStatistics(data, var = vars, grp = NULL, statistics = c("n", "mean", "sd", "min", "max"))
    #cat("\n")
  }
  
  if (doNormalTest) {
    normTestOut <- mshapiro.test(t(data[,vars]))
    #cat("MULTIVARIATE NORMALITY TEST\n")
    #cat("---------------------------\n")
    printTitle("MULTIVARIATE NORMALITY TEST")
    normTestTable <- as.table(cbind("Wilk-Shapiro", formatC(normTestOut$statistic[[1]],digits = 4, format="f"), formatC(normTestOut$p[[1]], digits = 4, format="f")))
    rownames(normTestTable) <- ""
    colnames(normTestTable) <- c("Statistic", "Value", "Prob")
    print(normTestTable, quote=FALSE)
  }
  
  if (doBoxM) {
    boxMTestOut <- BoxsMTest(data, vars, grpVar, alpha = 0.05)
    #cat("\nBOX'S M TEST FOR HOMOGENEITY OF COVARIANCE MATRICES\n")
    #cat("---------------------------------------------------\n")
    printTitle("BOX'S M TEST FOR HOMOGENEITY OF COVARIANCE MATRICES")
    
    if (!(is.nan(as.matrix(boxMTestOut$MBox)))) {
      if (length(boxMTestOut) == 4) {
        boxMTestTable<- as.table(cbind(formatC(boxMTestOut$MBox, digits = 4, format="f"), formatC(boxMTestOut$ChiSq, digits = 4, format="f"), formatC(boxMTestOut$df, digits = 0, format="f"), formatC(boxMTestOut$pValue, digits = 4, format="f")))
        rownames(boxMTestTable)   <- ""
        colnames(boxMTestTable) <- c("Box's M", "Chi-Square", "df", "Prob")
      } else {
        boxMTestTable<- as.table(cbind(formatC(boxMTestOut$MBox, digits = 4, format="f"), formatC(boxMTestOut$F, digits = 4, format="f"), formatC(boxMTestOut$df1, digits = 0, format="f"), formatC(boxMTestOut$df2, digits = 0, format="f"), formatC(boxMTestOut$pValue, digits = 4, format="f")))
        rownames(boxMTestTable)   <- ""
        colnames(boxMTestTable) <- c("Box M", "F", "df1", "df2", "Prob")
      }
      print(boxMTestTable)
    } else 
      cat("***\nBox's M statistic cannot be computed. \nDeterminant of pooled covariance matrix is less than 0.\n***")
    
  }
  
  cat("\n\nDISCRIMINANT ANALYSIS")

  discModel <- paste(grpVar, "~", paste(vars, collapse = " + "))
  
  if (is.null(priorProb)) {command <- paste("discOut <- lda(",discModel,", data = ", nameData,", const = TRUE)", sep = "")
  } else command <- paste("discOut <- lda(",discModel,", data = ", nameData,", prior = ", priorProb, ", const = TRUE)", sep = "")
  
  eval(parse(text = command))
 
  cat("\n\nSUBGROUP MEANS\n\n")
  print(round(discOut$means,4))
  
  #manova
  dataTest <- as.matrix(data[,vars])
  discOut.m <- manova(dataTest~data[,grpVar])
  discOut.sm <- summary(discOut.m,test="Wilks")
  
  cat("\n\nTEST OF SIGNIFICANCE OF THE DISCRIMINANT FUNCTION(S)\n\n")
  rownames(discOut.sm$stats)[1] <- grpVar
  print(discOut.sm)
   
  cat("\n\nDISCRIMINANT FUNCTION COEFFICIENTS\n\n")
  print(round(discOut$scaling, 4))
  cat("\n")
  
  #Classification Results
  if (is.null(priorProb)) { command2 <- paste("discOut2 <- lda(",discModel,", CV = TRUE, data = ",nameData,")", sep = "")
  } else command2 <- paste("discOut2 <- lda(",discModel,", CV = TRUE, data = ", nameData,", prior = ", priorProb, ")", sep = "")
  
  eval(parse(text = command2))
  #sort data by grpVar first
  
  newData <- cbind(data,discOut2$class)
  colnames(newData) <- c(names(data), "pClass")
  
  countsTable <- as.matrix(table(newData[,grpVar],newData[,"pClass"]))
  discOut$classResCountsTable <- cbind(countsTable, table(newData[,grpVar]))
  colnames(discOut$classResCountsTable) <- c(colnames(countsTable), "Total")

  cat("CLASSIFICATION RESULTS\n\nOriginal\\Predicted (Counts)\n")
  print(discOut$classResCountsTable)
  
  cat("\nOriginal\\Predicted (%)\n")
  pctsTable <- prop.table(countsTable,1)*100
  discOut$classResPctsTable <- cbind(pctsTable, rep(100, nrow(pctsTable)))
  colnames(discOut$classResPctsTable) <- c(colnames(pctsTable), "Total")
  print(round(discOut$classResPctsTable, 2))
  
  discOut$overallPct <- sum(countsTable[row(countsTable) == col(countsTable)])*100/sum(countsTable)
  cat("\n", round(discOut$overallPct, 2), "% of original cases correctly classified", sep = "")
  
  #classifying new observations
  if (classifyNew) {
    discOut3 <- predict(discOut, data2)
    data3 <- cbind(data2,discOut3$class)
    names(data3) <- c(names(data2), "Class")
    write.csv(data3, paste(outputPath,"discOut.csv", sep = ""))
    discOut$classified <- data3
  }

  discOut$signif <- discOut.sm
  
  return(invisible(list(Results = discOut)))
  
}  
  