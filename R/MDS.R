#---------------------------------------------------------#
#MDS - Function for performing Multidimensional Scaling   #
#
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# inputType - whether input is from raw data or distance Matrix
# vars - vector of names of the numeric variables
# idVar - name of the ID variable
# distance - distance measure to use
# type - type of MDS to perform
# dimnum - number of dimensions to be used
# useIdVar - logical; whether an ID variable will be used to label points in the plot or not
# pChars - vector representing characters for the points in the plot
# pSizes - vector representing sizes of points in the plot
# pCol - vector representing colors of points in the plot 
# showLeg - logical; whether a legend is displayed or not
# legTitle - title of the legend, if any
# legPos - position of the legend in the plot 
# legNcol - number of columns for the legend
#---------------------------------------------------------#

MDS <-function(data, outputPath, inputType = c("raw", "distMat"), vars = NULL, idVar, type = c("Classical", "Nonmetric"), 
               distClass = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Minkowski"),
               distNonmet = c("Manhattan", "Euclidean", "Canberra", "Bray", "Kulczynski", "Jaccard", "Gower", "AltGower", "Morisita", "Horn", "Mountford", "Raup", "Binomial", "Chao", "Cao"), 
               dimnum = 2, useIdVar = FALSE, pChars, pSizes, pCol, 
               showLeg = FALSE, legTitle, legPos = "bottomright", legNcol = 1, descriptive = FALSE, correlate = FALSE)
  
      UseMethod("MDS")

MDS.default <- function(data, outputPath, inputType = c("raw", "distMat"), vars = NULL, idVar, type = c("Classical", "Nonmetric"), 
                        distClass = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Minkowski"),
                        distNonmet = c("Manhattan", "Euclidean", "Canberra", "Bray", "Kulczynski", "Jaccard", "Gower", "AltGower", "Morisita", "Horn", "Mountford", "Raup", "Binomial", "Chao", "Cao"), 
                        dimnum = 2, useIdVar = FALSE, pChars, pSizes, pCol, 
                        showLeg = FALSE, legTitle, legPos = "bottomright", legNcol = 1, descriptive = FALSE, correlate = FALSE) 
  {
  
  if (is.character(data)) { 
		nameData <- data
	
    if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exist.", sep = "")) }
		tempData <- eval(parse(text = data)) 
  	
  } else {
      if (is.data.frame(data)) { 
		  nameData <- paste(deparse(substitute(data)))	
		  tempData <- data	
		  } else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
  }
  
  if (inputType == "raw") {
    data[, idVar] <- factor(data[, idVar])
    
    # -- PRINTING CLASS LEVEL INFORMATION -- #
    #ClassInformation(data[, c(vars, idVar)], respvar = vars)
    #cat("\n\n")
    
    numObsNM = nrow((na.omit(data[,vars])))
    numObs = nrow((data[,vars]))

    if (numObs == numObsNM) { cat("Number of Observations:", numObs,"\n")
    } else {
	    cat("Number of Observations:", numObs,"\n")
	    cat("Number of Observations Used:", numObsNM,"\n\n")
    }

    tempData <- na.omit(tempData[,c(vars, idVar)])  
    
  } else {
    tempData <- na.omit(tempData)
    d2 <- as.matrix(tempData[,-1])
#     rownames(d2) <- tempData[,1]
    rownames(d2) <- NULL
    colnames(d2) <- NULL
    if (!isSymmetric(d2)) { stop ("Input matrix should be symmetric.")}
  }
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
# 	if (!is.character(vars)) 	{ stop(paste("The object 'vars' should be a character vector.", sep = "")) }
# 	if (any(is.na(match(vars, names(tempData))))) { stop("At least one item in the character vector 'vars' does not match any variable name in the dataset.") }
    
  type <- match.arg(type)
	options(width = 5000, digits = 6)
#   distClass <- match.arg(distClass)
#   distNonmet <- match.arg(distNonmet)
  
  if (descriptive) { 
    DescriptiveStatistics(data = tempData, var = vars, statistics = c("nnmiss", "mean", "sd", "se.mean"))
    cat("\n\n")
  }	
  
  if (correlate) {
    BivariateCorrelationTest(data = tempData, var = vars, method = "pearson", alternative = "two.sided", statistics = FALSE)
    cat("\n\n")
  }
  
	if (type == "Classical") {
	  cat("CLASSICAL MULTIDIMENSIONAL SCALING \n\n")
    if (inputType == "raw") {
      distClass <- match.arg(distClass)
      d <- dist(tempData[,vars], method = tolower(distClass))
    } else d <- tempData[,-1]

    capture.output(fit <- cmdscale(d, eig=TRUE, k=as.numeric(dimnum))) # k is the number of dimensions

    points <- data.frame(fit$points)
    colnames(points) <- make.unique(rep("MDS",dimnum+1),sep ="")[2:(dimnum+1)]
    cat("\nPOINTS\n\n")
  	print(points)
		cat("\n")
		
	  eigenOut <- data.frame(fit$eig)
# 	  cat("\nEIGENVALUES\n")
# 	  colnames(eigenOut) <- ""
#     print(eigenOut)
    if (nrow(eigenOut)>10) { xLength = 10
    } else xLength = nrow(eigenOut)
	  png(filename = paste(outputPath,"MDS_Screeplot.png",sep=""))
	  plot(as.numeric(rownames(eigenOut))[1:xLength], eigenOut[1:xLength,1], xaxt = "n", xlab = "Eigenvalue number", ylab = "Eigenvalue", main = paste(type,"Scree Plot"), 
         type = "b", pch = 16)
    axis(1, at = as.numeric(rownames(eigenOut))[1:xLength])
	  dev.off()

    cat("\nP_", dimnum, " criterion: ", round(fit$GOF[2], digits = 4), sep = "") #2nd GOF statistic
	  
  	}else if (type == "Nonmetric"){
  	  cat("NON-METRIC MULTIDIMENSIONAL SCALING \n\n")
#   	  distNonmet <- match.arg(distNonmet)
  	  if (inputType == "raw") {
  	    distNonmet <- match.arg(distNonmet)
        if (tolower(distNonmet) == "bray") {
          d <- vegdist(tempData[,vars], method = "bray")
        } else d <- dist(tempData[,vars], method = tolower(distNonmet))
  	    capture.output(fit <- monoMDS(dist = d, k=as.numeric(dimnum)))
#       capture.output(fit <- metaMDS(tempData[, vars], distance = tolower(distNonmet), k=as.numeric(dimnum), model="global"))
  	  } else {
          d <- as.dist(tempData[,-1])
  	      capture.output(fit <- monoMDS(dist = d, k=as.numeric(dimnum)))
##        capture.output(fit <- metaMDS(tempData[,-1], distance = NULL, k=as.numeric(dimnum), model="global", autotransform = FALSE, wascores = FALSE, noshare = FALSE))
#       capture.output(fit <- metaMDS(tempData2, distance = distNonmetric, k=as.numeric(dimnum), model="global"))
  	  }
  		cat("\nPOINTS \n\n")
		  print(fit$points[1:dim(fit$points)[1], 1:dimnum])
		
#       #print gof
#   	  cat("\nGOODNESS OF FIT \n ")
#   	  gofOut <- data.frame(goodness(fit))
#       colnames(gofOut) <- ""
#       print(gofOut[1:length(gofOut)])

      #print stress statistic
      stressStat <- as.matrix(round(fit$stress, digits = 4))
  	  rownames(stressStat) <- ""
  	  cat("\nStress Statistic: ", stressStat[1,1], sep ="")
  	  cat("\n")
      
#       #create Shepard diagram
#       png(filename = paste(outputPath,"ShepardDiagram.png",sep=""))
#       shepOut <- Shepard(d=d,x=fit$points)
#   	  plot(shepOut, pch = ".")
#   	  lines(shepOut$x, shepOut$yf, type = "S")
#       dev.off()
  	  
      #create stressplot
  	  png(filename = paste(outputPath,"Stressplot.png",sep=""))
  	  stressplot(fit)
  	  dev.off()
  	        
#       #create scree plot
#       fit2<- NULL
#       stressValues <- NULL
#   	  if (nrow(tempData)>10) { nStress = 10
#   	  } else nStress = nrow(tempData)
# #       nStress = nrow(tempData)
#       for (i in 1:(nStress)) {
#         fit2 <- monoMDS(d=d, k=i)
#         stressValues[i] <- fit2$stress/100
#       }
      
#   	  png(filename = paste(outputPath,"MDS_Screeplot.png",sep=""))
#   	  plot(c(1:nStress), stressValues, xaxt = "n", xlab = "Dimension number", ylab = "Stress statistic", main = paste(type,"Scree Plot"), 
#   	       type = "b", pch = 16)
#       axis(1, at = c(1:nStress))
#   	  dev.off()
	}
  
  #create plot
  if (inputType == "raw")  {
    grpLevels = levels(factor(data[, idVar]))
    idData = data[,idVar]
  } else {
    grpLevels = levels(factor(data[,1]))
    idData = data[,1]
  }

#  if(!useIdVar) {
  grpCode<- data.frame(grpLevels, pChars, pCol)
  pCharNames <- make.unique(c(names(data),"pCharCode"), sep = "")
  pCharName <- pCharNames[length(pCharNames)]
  pColName <- make.unique(c(pCharNames,"pColCode"), sep = "")
  pColName <- pColName[length(pColName)]
  tempData2 <- cbind(data,grpCode$pChars[match(idData, grpCode[,"grpLevels"])], grpCode$pCol[match(idData, grpCode[,"grpLevels"])])
  colnames(tempData2)[ncol(tempData2)-1]<- pCharName
  colnames(tempData2)[ncol(tempData2)]<- pColName
#  } else {
#    pCharName = NULL
#    pColName = NULL
#    tempData2 = data
#  }
  mdsGraph(tempData2, idVar, outputPath, pCharName, pColName, grpLevels, idData, dimnum, fit, type, useIdVar, pChars, pSizes, pCol,
           showLeg, legTitle, legPos, legNcol)
  
	return(list(Results = fit, Method= type))
}#-- end stmt (mds Function)--#


