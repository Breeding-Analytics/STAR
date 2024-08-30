# --------------------------------------------------------#       
# MAOVTest: performs multivariate anova
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where output(s) will be saved
# yVars - vector of the numeric variables of interest
# factorVar - treatment factor
# repVar - blocking factor
# testStat - multivariate test statistic to use
# descriptive - logical; whether descriptive statistics will be displayed or not
# correlate - logical; whether pairwise correlations of yVars will be displayed or not
# doNormalTest - logical; whether multivariate normality test will be performed or not
# doBoxM - logical; whether Box's M test will be performed or not
# doSphericity - logical; whether Bartlett's Sphericity test will be performed or not
# saveSSH - logical; whether sums of squares and products for the hypothesis will be saved to a file or not
# saveSSE - logical; whether sums of squares and products for error will be saved to a file or not
# --------------------------------------------------------#

MAOVTest <- function(data, outputPath, yVars, factorVar, repVar = NULL, testStat = c("Pillai","Wilks","Hotelling-Lawley","Roy"), 
                     mymodel, descriptive = FALSE, correlate = FALSE, doNormalTest = FALSE, doBoxM= FALSE, doSphericity = FALSE, 
                     saveSSH = FALSE, saveSSE = FALSE)
  UseMethod("MAOVTest")

MAOVTest.default <- function(data, outputPath, yVars, factorVar, repVar = NULL, testStat = c("Pillai","Wilks","Hotelling-Lawley","Roy"), 
                     mymodel, descriptive = FALSE, correlate = FALSE, doNormalTest = FALSE, doBoxM= FALSE, doSphericity = FALSE, 
                     saveSSH = FALSE, saveSSE = FALSE) {

  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  tempData <- data
  
  #converts to factor the grouping variable(s)
  if (!is.null(repVar)) { tempData[,repVar] <- factor(tempData[,repVar]) }  
  if (!is.null(factorVar)) { tempData[,factorVar] <- factor(tempData[,factorVar]) }  
  
  # -- PRINTING CLASS LEVEL INFORMATION -- #
  ClassInformation(tempData[, c(factorVar, repVar, yVars)], respvar = c(factorVar, repVar, yVars))
  cat("\n\n")
  
  data <- na.omit(data[,c(yVars,factorVar,repVar)]) 
  
  #check if data is not empty
  if (nrow(data) == 0) { stop("Data set does not have any observation.") }
  
  #converts to factor the grouping variable(s)
  if (!is.null(factorVar)) { data[,factorVar] <- factor(data[,factorVar]) }  
  
  if (!is.null(repVar)) { 
    data[,repVar] <- factor(data[,repVar]) 
    nBalObs <- nlevels(data[,repVar]) * nlevels(data[,factorVar])
    if (nrow(data) != nBalObs) {stop("Data set should be balanced.")}
  }  
    
  results <- list()

  ###OPTIONS
  #DESCRIPTIVE     
  if (descriptive) { 
    DescriptiveStatistics(data, var = yVars, grp = NULL, statistics = c("n", "mean", "sd", "min", "max"))
    cat("\n")
  }

  #PAIRWISE CORRELATION 
  if (correlate) {
    cat("\n")
    BivariateCorrelationTest(data, var = yVars, method = "pearson", alternative = "two.sided", statistics = FALSE)
  }
  
  if (doNormalTest) {
    normTestOut <- mshapiro.test(t(data[,yVars]))
    #cat("MULTIVARIATE NORMALITY TEST\n") # hide by AAGulles 05.06.2015
    #cat("---------------------------\n") # hide by AAGulles 05.06.2015
    printTitle("MULTIVARIATE NORMALITY TEST") # added by AAGulles 05.06.2015
    normTestTable <- as.table(cbind("Wilk-Shapiro", formatC(normTestOut$statistic[[1]],digits = 4, format="f"), formatC(normTestOut$p[[1]], digits = 4, format="f")))
    rownames(normTestTable) <- ""
    colnames(normTestTable) <- c("Statistic", "Value", "Prob")
    print(normTestTable, quote=FALSE)
  }
  
  if (doBoxM) {
    boxMTestOut <- BoxsMTest(data, yVars, factorVar, alpha = 0.05)
    #cat("\n\nBOX'S M TEST FOR HOMOGENEITY OF COVARIANCE MATRICES\n") # hide by AAGulles 05.06.2015
    #cat("---------------------------------------------------\n") # hide by AAGulles 05.06.2015
    printTitle("BOX'S M TEST FOR HOMOGENEITY OF COVARIANCE MATRICES") # added by AAGulles 05.06.2015
    
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
  
  if (doSphericity)  {
    R <- cor(data[, yVars])
    n <- nrow(data)
    p <- length(yVars)
    chiSqVal <- -(n-1-(2*p+5)/6)*log(det(R))
    chiSqDf <- p*(p-1)/2
    pValue = pchisq(chiSqVal, chiSqDf, lower.tail = F)
    
    #cat("\n\nBARTLETT'S SPHERICITY TEST\n") # hide by AAGulles 05.06.2015
    #cat("--------------------------\n") # hide by AAGulles 05.06.2015
    printTitle("BARTLETT'S SPHERICITY TEST") # added by AAGulles 05.06.2015
    sphereTestTable <- as.table(cbind("Chi-Square", formatC(chiSqVal, digits = 4, format="f"), formatC(chiSqDf, digits = 0, format="f"), formatC(pValue, digits = 4, format="f")))
    rownames(sphereTestTable) <- ""
    colnames(sphereTestTable) <- c("Statistic", "Value", "Df", "Prob")
    print(sphereTestTable, quote=FALSE)
  }
  ###end - OPTIONS
  
#   modelLHS = paste("cbind(", paste(yVars, collapse=", "), ")", sep = "")
#   modelRHS = paste(factorVar, sep = "")
#   if (!is.null(repVar)) modelRHS = paste(factorVar, repVar, sep = " + ")
#   mymodel = paste(modelLHS, " ~ ", modelRHS, sep = "")
  linmodel = (lm((mymodel), data = data))

  maov  <- manova(linmodel, data = data)

  #cat("\n\nSUM OF SQUARES AND CROSS PRODUCTS FOR THE HYPOTHESIS \n") # hide by AAGulles 05.06.2015
  #cat("----------------------------------------------\n") # hide by AAGulles 05.06.2015
  printTitle("SUM OF SQUARES AND CROSS PRODUCTS FOR THE HYPOTHESIS") # added by AAGulles 05.06.2015

  cat(paste("\nTERM:", factorVar, "\n\n"))
  trtSS = summary(maov)$SS[[factorVar]]
  print(trtSS)
  if (saveSSH)
    write.csv(trtSS, paste(outputPath,"/", factorVar, "SS.csv", sep = ""))
  if (!is.null(repVar)) {
    cat(paste("\nTERM:", repVar, "\n\n"))
    repSS = summary(maov)$SS[[repVar]]
    print(repSS)
    if (saveSSH)
      write.csv(repSS, paste(outputPath,"/", repVar, "SS.csv", sep = ""))
  }  
  
  #cat("\n\nSUM OF SQUARES AND CROSS PRODUCTS FOR ERROR\n") # hide by AAGulles 05.06.2015
  #cat("-------------------------------------\n") # hide by AAGulles 05.06.2015
  printTitle("SUM OF SQUARES AND CROSS PRODUCTS FOR ERROR") # added by AAGulles 05.06.2015
  
  errSS = summary(maov)$SS$Residuals
  print(errSS)
  if (saveSSE)
    write.csv(errSS, paste(outputPath,"/errorSS.csv", sep = ""))
  
  if (is.null(repVar)) { design = "CRD"
  } else design = "RCBD"
  
  #cat(paste("\n\nMULTIVARIATE ANALYSIS OF VARIANCE IN ", design, "\n", sep = "")) # hide by AAGulles 05.06.2015
  #cat("-----------------------------------------\n") # hide by AAGulles 05.06.2015
  printTitle(paste("MULTIVARIATE ANALYSIS OF VARIANCE IN ", design, sep = "")) # added by AAGulles 05.06.2015

  print(summary(maov, test = testStat))
                      
}