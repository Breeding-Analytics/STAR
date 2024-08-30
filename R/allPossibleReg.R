#----------------------------------------------------#
# allPossibleReg - function for performing all possible regression
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

allPossibleReg <- function(data, depVar, indepVar, selectionStat = "mallow", constant = TRUE) 
  UseMethod("allPossibleReg")

allPossibleReg.default <- function(data, depVar, indepVar, selectionStat = "mallow", constant = TRUE) {
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  for (m in 1:length(depVar)) {
    #missing values?
    data <- na.omit(data)
    
    ncombn <- vector()
    nIndepVar = length(indepVar)
    for (k in 1:nIndepVar) {
      ncombn[[k]] <- dim(combn(nIndepVar,k))[2]
    }
    bigN <- max(ncombn)
    regModel <- paste(depVar[m] , "~", paste(indepVar, collapse = " + ")) 
    regCommand <- paste("reg1 <- regsubsets(", regModel, ", data = data, nbest = ", bigN, ", intercept = ", constant, ", really.big = TRUE)", sep = "")
    eval(parse(text = regCommand))
    sumReg1<-summary(reg1)
    selStat <- function(meas) {
      switch(meas, 
           mallow = sumReg1$cp,
           rsquare = sumReg1$rsq,
           adjrsq = sumReg1$adjr2)
    }
    selStatName <- function(meas) {
      switch(meas, 
             mallow = "Mallow's C(p)",
             rsquare = "R-Square",
             adjrsq = "Adjusted R-Square")
    }
    selStats <- selStat(selectionStat)
    startInd = 1
    desStats <- vector()
    desInd <- vector()
    selIndepVar <- vector()
    regResult <- NULL
    for (o in 1:nIndepVar){ 
      batchEnd <- startInd+ncombn[o]-1
      #determine best "model" based on statistic
      if (selectionStat == "mallow") {
        desStats[[o]] = min(abs(selStats[startInd:batchEnd] - o))
        desInd[[o]] <- which(abs(selStats[startInd:batchEnd] - o) == desStats[[o]]) + startInd - 1
      } else { 
        desStats[[o]] <- max(selStats[startInd:batchEnd])
        desInd[[o]] <- which(selStats[startInd:batchEnd] == desStats[[o]]) + startInd - 1}
      startInd <- startInd + ncombn[o]
      if (constant) { selVarsInd <- which(sumReg1$which[desInd[[o]],-1])
      } else selVarsInd <- which(sumReg1$which[desInd[[o]],])
      
      selIndepVar[[o]] <- paste(names(selVarsInd),collapse=" + ")
      dataO <- data[,c(depVar[m], names(selVarsInd))]
      if (constant) { regCommand <- paste("regResult[[", o, "]] <- lm(", depVar[m], " ~ ", selIndepVar[[o]],", data = dataO)", sep = "")
      } else { regCommand <- paste("regResult[[", o, "]] <- lm(", depVar[m], " ~ ", selIndepVar[[o]], "-1, data = dataO)", sep = "")}
      
      eval(parse(text = regCommand))
    }
    selCp <- format(sumReg1$cp[desInd],digits=4, justify = "centre")
    selRSq <- format(sumReg1$rsq[desInd]*100,digits=4, justify = "centre") #in pct
    selAdjRSq <- format(sumReg1$adjr2[desInd]*100,digits=4, justify = "centre") #in pct
    selTable <- as.data.frame(cbind(format(selIndepVar,justify = "left"), selCp, selRSq, selAdjRSq))
    
    colnames(selTable) = c(format("Variables", width = nchar(max(selIndepVar)), justify = "centre"), 
                           format("C(p)", width = nchar(max(selCp)), justify = "centre"),
                           format("R-square", width = nchar(max(selRSq)), justify = "centre"),
                           format("Adj R-sq", width = nchar(max(selRSq)), justify = "centre"))
    
    cat("\nResponse Variable: ", depVar[m], "\n")
    cat("\nAll-possible regression (", selStatName(selectionStat), ")\n", sep = "")
    cat("\nSummary Table: \n")
    printDataFrame(selTable)
    cat("\n")
    
    for (o in 1:length(indepVar)) {
     cat("Model ", o, ": ", sep="")
     cat("\nIndependent Variables: ", selIndepVar[o])
     cat("\n\nParameter Estimates:\n", sep="")
     
     coefTable <- data.frame(rownames(summary(regResult[[o]])$coefficient),summary(regResult[[o]])$coefficient)
     colnames(coefTable) <- c("Variable", attributes(summary(regResult[[o]])$coefficient)$dimnames[[2]])
     coefTable[,1] <- as.character(coefTable[,1])
     if (constant) { coefTable[1,1] <- "Intercept" }
     printDataFrame(coefTable)
     cat("\n")
     
    }    
  }
}