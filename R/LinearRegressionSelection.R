#----------------------------------------------------#
# LinearRegressionSelection - function for performing Linear Regression Analysis with a selection method
#
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# depVar - vector of names of the response variable(s)
# indepVar - vector of names of the predictor variable(s)
# constant - logical; whether constant is included in the model or not
# correlate - logical; whether parametric correlation coefficients will be estimated for all pairs of variables
# selection - selection method to use, if any
# selectionStat - statistic to be used, if selection procedure used is all possible regression
# confInt - logical; whether confidence interval estimates will be constructed for the coefficients or not
# confLevel - level of confidence coefficient, if confidence interval is to be constructed
# covMatrix - logical; whether the variance-covariance matrix will be displayed or not
# normality - vector of names of tests for normality to be performed, if any
# heteroskedasticity - vector of names of tests for heteroskedasticity to be performed, if any
# autoCorr - logical; whether test for autocorrelation will be performed or not
# VIF - logical; whether variance inflation factor will be computed or not
# COOKS - logical; whether cooks-weisberg test will be performed or not
# leverage - logical; whether leverage values will be computed or not
#----------------------------------------------------#

LinearRegressionSelection <- function(data, depVar, indepVar, constant = TRUE, correlate = FALSE, 
                                     selection = "none", selectionStat = NULL,
                                     confInt = FALSE, confLevel = 0.95, 
                                     covMatrix = FALSE, normality = NULL, 
                                     heteroskedasticity = NULL, autoCorr = FALSE, 
                                     VIF = FALSE, COOKS = FALSE, leverage = FALSE) 
  UseMethod("LinearRegressionSelection")
  
LinearRegressionSelection.default <-function(data, depVar, indepVar, constant = TRUE, correlate = FALSE, 
                                             selection = "none", selectionStat = NULL,
                                             confInt = FALSE, confLevel = 0.95, 
                                             covMatrix = FALSE, normality = NULL, 
                                             heteroskedasticity = NULL, autoCorr = FALSE, 
                                             VIF = FALSE, COOKS = FALSE, leverage = FALSE) {
  
  if (is.character(data)) {
    nameData <- data
    data <- eval(parse(text = data))
  } else { nameData <- paste(deparse(substitute(data))) } #reads data: if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #descriptive statistics
  DescriptiveStatistics(data = data, var = c(depVar,indepVar), statistics = c("nnmiss", "mean", "sd", "se.mean")) 
  cat("\n")

  #correlation
  if (correlate) {
    cat("\n")
    BivariateCorrelationTest(data = data, var = c(depVar, indepVar), method = "pearson", alternative = "two.sided", statistics = FALSE)
  }
  
  regOut = NULL
  outData = NULL
  outFit = NULL
  
  #selection
  if (selection == "none") {
    #perform regression
    regOut = LinearRegressionAnalysis(data, depVar, indepVar, constant, statistics = FALSE, confInt, confLevel, covMatrix,
                             normality, heteroskedasticity, autoCorr, VIF, COOKS, leverage)
    
  } else if (selection == "allposs") {
    allPossibleReg(data, depVar, indepVar, selectionStat, constant)
  
  } else {
    selMethod <- switch(selection, 
           forward = "Forward selection",
           backward = "Backward elimination",
           stepwise = "Stepwise regression")
    
    
    cat("SELECTION METHOD: ", selMethod, "\n\n", sep = "")
      if (selection == "stepwise") selection = "both"
    
    for (i in 1:length(depVar)) {
      cat("DEPENDENT VARIABLE: ", depVar[i], "\n\n", sep = "")
      
      if (constant) {regCom <- paste("lm(", depVar[i], " ~ ", paste(indepVar, collapse = "+"),", data = data)", sep = "")
      } else { regCom <- paste("lm(", depVar[i], " ~ ", paste(indepVar, collapse = "+"),"- 1, data = data)", sep = "")}
      res <- eval(parse(text = regCom))
      fwdReg <- step(res, direction = selection, trace = 1)
        
      if (constant) { selIndepVar <- names(fwdReg$coefficients)[-1]
      } else selIndepVar <- names(fwdReg$coefficients)
      
      cat("\n")
      #perform regression
      regOut = LinearRegressionAnalysis(data, depVar[i], selIndepVar, constant, statistics = FALSE, confInt, confLevel, covMatrix,
                                 normality, heteroskedasticity, autoCorr, VIF, COOKS, leverage)
      if (i == 1) { outData = regOut$data
      } else {
        outVarNames = NULL
        outVarNames[[1]] = paste("regOut$data$", depVar[i], "_pred",sep = "")
        outVarNames[[2]] = paste("regOut$data$", depVar[i], "_resid", sep = "")
        outVars = data.frame(cbind(eval(parse(text = outVarNames[[1]])), eval(parse(text = outVarNames[[2]]))))
        colnames(outVars) = c(paste(depVar[i], "_pred",sep = ""), paste(depVar[i], "_resid",sep = ""))
        outData = cbind(outData, outVars)
      }
      outFit[i] = regOut$modelFit
      
    }      
      
  }
        

  if (selection == "allposs") {
    return(list(data = data))
  } else if (selection == "none") {
    return(list(data = regOut$data, modelFit = regOut$modelFit))
  } else {
    return(list(data=outData, modelFit = outFit))
  }

}
  
