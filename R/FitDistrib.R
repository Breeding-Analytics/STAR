# --------------------------------------------------------       
# FitDistrib: fits selected distribution to data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where output(s) will be saved
# yVar - numeric variable of interest
# byVar
# distrib - distribution to fit
# fitMethod - method used for fitting
# --------------------------------------------------------

FitDistrib <- function(data, outputPath, yVar, byVar = NULL, 
                       distrib = c("beta", "exp", "gamma", "geom", "lnorm", "logis", "nbinom", "norm", "pois"), 
                       fitMethod = c("mle", "mme", "mge"))
  UseMethod("FitDistrib")

FitDistrib.default <- function(data, outputPath, yVar, byVar = NULL, 
                       distrib = c("beta", "exp", "gamma", "geom", "lnorm", "logis", "nbinom", "norm", "pois"), 
                       fitMethod = c("mle", "mme", "mge")) {
 
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #creates a grouping variable with ones for all rows if no grouping variable is declared
  if (is.null(byVar)) { 
    grp.Var = rep(1,nrow(data))
    byVar = "grp.Var"
    data = cbind(data,grp.Var)
  }  
  
  #converts to factor the grouping variable(s)
  if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }  
  
  results = list()

  switch(fitMethod,
         mle={fit.Method = "maximum likelihood estimation"},
         mme={fit.Method = "moment matching estimation"},
         mge={fit.Method = "maximum goodness-of-fit estimation"},
         stop("Invalid fitting method.")
  )
  
  cat("DISTRIBUTION FITTING\nMETHOD: ", fit.Method, "\n\nVARIABLE: ", yVar, "\n\n", sep="")
  
  if (is.null(byVar)) numRun = 1  else numRun = nlevels(data[,byVar])
  
  for (j in 1:numRun) {
    results[[j]] <- list()
    
    #creates data by subgroup, if a grouping variable is defined
    if (!is.null(byVar)) {
      byVarLevel = levels(data[,byVar])[j]
      tempData = data[which(data[,byVar] == byVarLevel),]
      if (length(byVarLevel) > 1)  {
        cat(byVar, " = \"", byVarLevel, "\"\n", sep="")
        subTitle = paste(byVar,"=",byVarLevel, sep=" ")
      }
        
    }
    
    tempData1 <- tempData[(is.na(tempData[, yVar]) == FALSE),]
    tempData <- tempData1[,yVar]

    if (length(tempData) > 0) {
      distribn = NULL 
      
      allDiscrete = TRUE
      #fits the specified distribution(s)  
      for(i in 1:length(distrib)) {
        results[[j]][[i]] <- list()
        
        switch(distrib[i],
               beta = {distribn[i] = "Beta"},
               exp = {distribn[i] = "Exponential"},
               gamma = {distribn[i] = "Gamma"}, 
               geom  = {distribn[i] = "Geometric"}, 
               lnorm = {distribn[i] = "Log-normal"},
               logis = {distribn[i] = "Logistic"},
               nbinom = {distribn[i] = "Negative binomial"},
               norm = {distribn[i] = "Normal"},
               pois = {distribn[i] = "Poisson"},
               stop("Invalid distribution selected.")
        )
        
        results$distrib[[i]] = distribn[i]
        
        #check if all distributions are discrete
        if (is.element(distrib[i], c("beta", "exp", "gamma", "lnorm", "logis", "norm"))) allDiscrete <- FALSE
        
        cat("\nDISTRIBUTION: ", distribn[i], "\n\n", sep="")
        if (fitMethod == "mge") {
          results[[j]][[i]]$estimates = fitdist(data = tempData, distr = distrib[i], method = "mge", gof = "CvM")
        } else results[[j]][[i]]$estimates = fitdist(data = tempData, distr = distrib[i], method = fitMethod)
        print(results[[j]][[i]]$estimates)
        
        #goodness of fit test
        capture.output(results[[j]][[i]]$fitComp <- gofstat(results[[j]][[i]]$estimates))
        
        if (is.element(distrib[i], c("beta","exp", "gamma", "lnorm", "logis", "norm")))  {
          cat("\nGoodness-of-fit\n")
          gofComp = format(rbind(results[[j]][[i]]$fitComp$cvm, results[[j]][[i]]$fitComp$ad, results[[j]][[i]]$fitComp$ks), digits=4)
          gofDecisions = rbind(
            ifelse(!is.null(results[[j]][[i]]$fitComp$cvmtest), results[[j]][[i]]$fitComp$cvmtest, "NULL"), 
            ifelse(!is.null(results[[j]][[i]]$fitComp$adtest), results[[j]][[i]]$fitComp$adtest, "NULL"), 
            ifelse(!is.null(results[[j]][[i]]$fitComp$kstest), results[[j]][[i]]$fitComp$kstest, "NULL"))
          gofStats = cbind(gofComp, gofDecisions) 
          
          colnames(gofStats) <- c("Statistic", "Decision")
          rownames(gofStats) <- c("Cramer-von Mises", "Anderson-Darling", "Kolmogorov-Smirnov")
          results[[j]][[i]]$gofTable <- as.table(gofStats)
        } else {
          #chisquare
          cat("\nChi-squared Goodness-of-fit\n")
          gofChisq = cbind(
            ifelse(!is.null(results[[j]][[i]]$fitComp$chisq), format(results[[j]][[i]]$fitComp$chisq, digits = 4), "NULL"), 
            ifelse(!is.null(results[[j]][[i]]$fitComp$chisqpvalue), format(results[[j]][[i]]$fitComp$chisqpvalue, digits = 4), "NULL"), 
            ifelse(!is.null(results[[j]][[i]]$fitComp$chisqdf), results[[j]][[i]]$fitComp$chisqdf, "NULL")
          )
          colnames(gofChisq) <- c("Statistic", "p-value", "df")
          rownames(gofChisq) <- ""
          results[[j]][[i]]$gofTable <- as.table(gofChisq)
        }
        
        print(results[[j]][[i]]$gofTable)
        
      }
      
      #plots
      if (numRun > 1) { plotfile = paste(outputPath,"plots_", byVarLevel, ".png",sep="")     
      } else { plotfile = paste(outputPath,"plots.png",sep="")}
      
      png(filename = plotfile)
      
      estimatesList = list()
      for (i in 1:length(distrib))
        estimatesList[[i]] = results[[j]][[i]]$estimates
      
      if (allDiscrete) {
        par(mfrow=c(1, 2))
       plotDiscDist(tempData, distrib, estimatesList)   
        
      } else {
        if (length(distrib) > 1) {
          par(mfrow=c(2, 2))
          histData =  hist(tempData, plot = F)  
          minDensity = min(histData$density, na.rm = TRUE) 
          maxDensity = max(histData$density, na.rm = TRUE)

          cdfcomp(estimatesList, xlim = c(min(tempData), max(tempData)), ylim = c(0,1), legendtext=c(distrib))
          denscomp(estimatesList, ylim = c(minDensity, maxDensity), xlim = c(min(tempData), max(tempData) * 1.25), legendtext=c(distrib)) #ylim = c(0,0.0004),
          qqcomp(estimatesList, xlim = c(min(tempData)-0.25*min(tempData), max(tempData) * 1.25), legendtext=c(distrib))
          ppcomp(estimatesList, xlim = c(0,1), ylim = c(0,1), legendtext=c(distrib))
        } else {
           plotdist(tempData, distrib[1], c(as.list(estimatesList[[1]]$estimate))) 
        }
      }
      dev.off()
      cat("\n----------------------------------------\n\n")

      numDiscrete <- sum(is.element(distrib, c("geom", "nbinom", "pois")) == TRUE)
      if (length(distrib) > 1 && (numDiscrete > 0 && numDiscrete < length(distrib))) 
        cat("\nWarning: Variables are considered continuous in this function\n\n")
    } else {# end of if (length(tempData) > 0) 
      cat("\nError: No valid observations.\n\n")
    }
    
  }
  results$method = fit.Method
    
  return(results)
}