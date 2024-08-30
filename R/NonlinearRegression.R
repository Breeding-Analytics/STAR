######
# NonLinearRegression - function for performing NonLinear Regression Analysis
#
# data - name of R dataframe
# depVar - name of response variable
# indepVar - name of predictor variable
# statistics - logical; whether descriptive statistics are displayed for all variables or not
# model - type of model to fit
# modelExp - expression for the right-hand side of the model, if customized
# startValues - starting values for parameters of the model, if customized
# covMatrix - logical; whether variance-covariance matrix for coefficients will be printed or not
# normality - logical; whether wilk-shapiro test for normality is to be done or not
# confInterval - logical; whether confidence intervals for the parameters will be constructed or not
# confLevel - level of confidence for confidence interval to be constructed
######

NonLinearRegression <- function(data, depVar, indepVar, statistics = FALSE, model, modelExp = NULL, startValues = NULL,
    normality = FALSE, covMatrix = FALSE, confInterval = FALSE, confLevel = 0.95
)  UseMethod("NonLinearRegression")   

NonLinearRegression.default <- function(data, depVar, indepVar, statistics = FALSE, model, modelExp = NULL, startValues = NULL,
    normality = FALSE, covMatrix = FALSE, confInterval = FALSE, confLevel = 0.95
   ) {

  if (is.character(data)) {
    nameData <- data
    data <- eval(parse(text = data))
  } else { nameData <- paste(deparse(substitute(data))) } 

  #generates descriptive statistics
  if (statistics) {
    DescriptiveStatistics(data = data, var = c(depVar, indepVar), statistics = c("nnmiss", "mean", "sd", "se.mean")) 
    cat("\n")
  }  
  
  #defines model right-hand side
  getModelRHS <- function(iVar, model) {
    switch(model,
           quadratic = paste(iVar, "+ I(", iVar, "^2)", sep = ""),
           cubic = paste(iVar, "+ I(", iVar, "^2)", "+ I(", iVar, "^3)", sep = ""),
           logis = paste("SSlogis(", iVar, ", Asym, xmid, scal)", sep = ""),
           fpl = paste("SSfpl(", iVar, ", A, B, xmid, scal)",sep = ""),
           gompertz = paste("SSgompertz(", iVar, ", Asym, b2, b3)",sep = ""),
           weibull = paste("SSweibull(", iVar, ", Asym, Drop, lrc, pwr)",sep = ""),
           custom = modelExp
    )
  }
  
  #defines command for nonlinear regression
  getNlRegCom <- function(modToFit, dataName, startVals = NULL) {
    switch(modGrp,
           first = paste("nls1 <- lm(", modToFit, ", data =  ", dataName, ")", sep = ""),
           second = paste("nls1 <- nls(", modToFit, ", data = ", dataName, ")", sep = ""),
           third = paste("nls1 <- nls(", modToFit, ", data = ", dataName, ", start = ", startVals, ")", sep = "")
    )
  }
  

  nlsOut = NULL
  sumNls = NULL
  confOut = NULL
  data2 = NULL
  
  cat("NONLINEAR REGRESSION ANALYSIS\n\n")
  
  for (i in 1:length(depVar)) {

    cat("Dependent Variable:", depVar[i], "\n\n")
    
    if (model == "quadratic" || model == "cubic") { modGrp = "first"
    } else if (model == "custom") { modGrp = "third"
    } else modGrp = "second"
    
    
    modelRHS <- getModelRHS(indepVar, model)
    startValuesL <- startValues
    modelToFit<- paste(depVar[i], " ~ ", modelRHS)
    nlRegCom <- getNlRegCom(modelToFit, nameData, startValues)
    res <- try(eval(parse(text = nlRegCom)), silent = TRUE)
    if (class(res) == "try-error") {
      msg <- trimStrings(strsplit(res, ":")[[1]])
      msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
      msg <- gsub("\"", "", msg)
      cat("\nError in NonLinearRegressionAnalysis:\n",msg, "\n\n", sep = "")
      next
    }
    nlsOut[[i]] = nls1
    sumNls[[i]] = summary(nls1)
    
    #performs Wilk-Shapiro Test for Normality
    if (normality) {
      residData <- data.frame(sumNls[[i]]$residuals)
      colnames(residData) <- "residual"
      NormalityTest(data = residData, var = "residual", grp = NULL, method = "swilk") #swilk only
      cat("\n")
    }
  
    cat(paste("Formula: ", modelToFit, "\n\n", sep = ""))

    cat("Parameter Estimates: \n")
    confOut[[i]] <- try(confint2(nlsOut[[i]],level = confLevel), silent = TRUE)
    
    #creates table for point estimates and (if any) confidence interval estimates
    if (confInterval && (class(confOut[[i]]) != "try-error")) {
      if (modGrp == "first") { rownames(sumNls[[i]]$coefficient)[1] <- "Intercept" }
      if (nrow(sumNls[[i]]$coefficient) != 1) {
        coefTable <- data.frame(rownames(sumNls[[i]]$coefficient),sumNls[[i]]$coefficient[,1:3], confOut[[i]], sumNls[[i]]$coefficient[,4])     
      } else {
        coefTable <- data.frame(rownames(sumNls[[i]]$coefficient),t(sumNls[[i]]$coefficient[,1:3]), confOut[[i]], sumNls[[i]]$coefficient[,4])     
      }
      
      colnames(coefTable) <- c("Variable", attributes(sumNls[[i]]$coefficient)$dimnames[[2]][1:3], "LL CI*","UL CI*", attributes(sumNls[[i]]$coefficient)$dimnames[[2]][4])
      printDataFrame(coefTable)
      cat("* At ",(confLevel)*100,"% Confidence Interval\n\n", sep = "")
    } else {
      coefTable <- data.frame(rownames(sumNls[[i]]$coefficient),sumNls[[i]]$coefficient)
      colnames(coefTable) <- c("Variable", attributes(sumNls[[i]]$coefficient)$dimnames[[2]])
      printDataFrame(coefTable)
      cat("\n")
      if (confInterval && (class(confOut[[i]]) == "try-error"))
       cat("Error in constructing confidence interval: \n\n", msg, sep = "")
    }
    
    #computes and prints residual sum of squares
    rSS = sumNls[[i]]$df[2]*sumNls[[i]]$sigma^2
    cat(paste("Residual Sum of Squares: "), format(rSS, nsmall = 4), "\n\n", sep = "")
    
    #displays variance-covariance matrix
    if (covMatrix) {
      cat("Coefficient Variance-Covariance Matrix\n")
      vcovTable = vcov(nlsOut[[i]])
      if (modGrp == "first") {
        rownames(vcovTable)[1] <- "Intercept"
        colnames(vcovTable)[1] <- "Intercept"
      }
      print(vcovTable, digits = 5)
      cat("\n\n")
    }
    
  if (!is.null(nlsOut[[i]])) {
    # saves fitted values and residuals
    if (i == 1) {
      newVar1 <- make.unique(c(names(data),paste(depVar[i],"pred", sep = "_")), sep = "")
      data2 <- cbind(data, fitted(nlsOut[[i]]))      
      colnames(data2) <- newVar1
    } else {
      newVar1 <- make.unique(c(names(data2),paste(depVar[i],"pred", sep = "_")), sep = "")
      data2 <- cbind(data2, fitted(nlsOut[[i]]))      
      colnames(data2) <- newVar1
    }
    
    newVar2 <- make.unique(c(names(data2),paste(depVar[i],"resid", sep = "_")), sep = "")
    data2 <- cbind(data2, residuals(nlsOut[[i]]))
    colnames(data2) <- newVar2
    
    }  
  } # end of for (i in 1:length(depVar))
                          
  if (!is.null(data2)) data = data2
  
  return(invisible(list(nlsFit = sumNls, confInt = confOut, data = data)))
}