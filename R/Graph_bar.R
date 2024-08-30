# --------------------------------------------------------
# GraphBar: creates a barchart of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# nVar - name of the numeric variables
# cVar - name of the categorical variable
# clustVars - name of  variable used for creating cluster of bars
# mTitle - main title for bar graph, if any
# yAxisLab - labels to be used for the y-axis (for the numeric variables)
# xAxisLab - label to be used for x-axis (for the categorical variable)
# minValue - minimum value for the the Value Axis
# maxValue - maximum value for the Value Axis
# typeData - whether dataframe contains raw data or summary statistics
# descStat - whether the bars to be created will represent means, medians, sums,
#             or frequencies (if typeData is raw)
# barsHoriz - logical; whether bars are drawn horizontally or not
# barsClust - logical; whether bars from different response variables are clustered or not (stacked)
# byVar - name of categorical variable for which separate graphs are created
# errBars - logical; whether error bars are displayed or not
# typeErrBar - whether error bars displayed (if any) are based on:
# 					standard deviation, standard error, or confidence interval
# errMult - multiplier to be used on standard deviaton/standard error 
# confLev - confidence level, if confidence interval type error bars are to be displayed
# axisLabelStyle - style for the axes labels
# barColor - vector of rgb values for the color(s) of the bars (or names of colors)
# showLeg - logical; whether the graph legend is displayed or not
# legPos - position of the legend (if showLeg is TRUE)
# legTitle - text for the title of the legend, if displayed
# legCol - number of columns used for legend text
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# showCatVarLevels - logical; whether levels of the categorical variable will be displayed or not
# barDensity - number of lines used for shading the bars
# barLineAngle - vector of angle values for the lines inside the bars
# --------------------------------------------------------

GraphBar <- function(data, outputPath, nVar, cVar, clustVars = NULL, mTitle = NULL, yAxisLab = NULL, 
                     xAxisLab = NULL, minValue = NULL, maxValue = NULL, typeData = c("raw", "sumstat"), 
                     descStat = c(NULL, "mean", "median", "sum", "freq"), barsHoriz = FALSE, barsClust = TRUE,
                     byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdDev", "stdErr", "confInt"), 
                     errMult = 1, confLev = NULL, axisLabelStyle = 1, barColor = NULL, 
                     showLeg = FALSE, legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",
                                                 "center"), legTitle = NULL, legCol = 1, boxed = TRUE,
                     multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1,
                     orientGraphs = c("left-right", "top-bottom"),
                     showCatVarLevels = TRUE, barDensity, barLineAngle) 
  UseMethod("GraphBar")

GraphBar.default <- function(data, outputPath, nVar, cVar, clustVars = NULL, mTitle = NULL, yAxisLab = NULL, 
								xAxisLab = NULL, minValue = NULL, maxValue = NULL, typeData = c("raw", "sumstat"), 
								descStat = c(NULL, "mean", "median", "sum", "freq"), barsHoriz = FALSE, barsClust = TRUE,
								byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdDev", "stdErr", "confInt"), 
								errMult = 1, confLev = NULL, axisLabelStyle = 1, barColor = NULL, 
								showLeg = FALSE, legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",
										"center"), legTitle = NULL, legCol = 1, boxed = TRUE,
								multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1,
								orientGraphs = c("left-right", "top-bottom"),
								showCatVarLevels = TRUE, barDensity, barLineAngle) { 
	
	#reads data
	if (is.character(data)) { data <- eval(parse(text = data)) }
	
	#creates a grouping variable (clustVars) if it is not declared
  clustVar = NULL
	if (is.null(clustVars)) { 
	  clustVar <- make.unique(c(names(data), "clustVar"))[length(make.unique(c(names(data), "clustVar")))]
	  data[,clustVar] = rep(1,nrow(data))
	} else {
	  #combines levels of 2 or more line variables
	  if (length(clustVars) == 2) { 
	    clustVar = paste(clustVars[1], clustVars[2], sep = "-")
	    data[,clustVar] = paste(data[,clustVars[1]], data[,clustVars[2]], sep ="-")
	  } else if (length(clustVars) == 3) { 
	    clustVar = paste(clustVars[1], clustVars[2], clustVars[3], sep = "-")
	    data[,clustVar] = paste(data[,clustVars[1]], data[,clustVars[2]], data[,clustVars[3]], sep ="-")
	  } else {
      clustVar = clustVars[1]
      data[,clustVar] = data[,clustVars[1]]
	  }  
	}
	
	#converts to factor the grouping variable(s)
	if (!is.null(cVar)) { data[,cVar] <- factor(data[,cVar]) }
	if (!is.null(clustVar)) { data[,clustVar] <- factor(data[,clustVar]) }
  if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }	
	
	if (!multGraphs) {
		numRowsGraphs = 1
		numColsGraphs = 1
	} 
	
	#determines number of cells allocated for graphs (esp. for multiple graphs)
 	numCells = numRowsGraphs * numColsGraphs
	
	numGroups = 1
  
	#determines number of graphs to be created
	if (!is.null(byVar)) {
	  numGroups = nlevels(data[,byVar])
	  numGraphs = nlevels(data[,byVar])*length(nVar)
	} else numGraphs = length(nVar)

 	graphNum = 1
	#counts the number of files to save
	k = 1
	
	for (m in 1:numGroups) { 
    
	  #creates data by subgroup, if a grouping variable is defined
	  if (!is.null(byVar)) {
	    tempData1 = data[which(data[,byVar] == levels(data[,byVar])[m]),]
	    subTitle = paste(byVar,"=",levels(tempData1[,byVar])[m], sep=" ")
	  } else {
	    tempData1 = data
	    subTitle = NULL
	  }
	  
	  for (i in 1:length(nVar)) {
	    #creates device for saving graph(s)
	    if (!multGraphs) {
	      png(filename = paste(outputPath,"bargraph_",k,".png",sep=""))
	      par(mfrow=c(numRowsGraphs, numColsGraphs))
	      
	    } else {
	      if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)	{
          widthAdj = numColsGraphs * 480
          heightAdj = numRowsGraphs * 480
	        png(filename = paste(outputPath,"bargraph_",k,".png",sep=""), width = widthAdj, height = heightAdj)
	        
	        if (orientGraphs == "top-bottom") {
	          par(mfcol=c(numRowsGraphs, numColsGraphs))
	        } else if (orientGraphs == "left-right") {
	          par(mfrow=c(numRowsGraphs, numColsGraphs))
	        }
	      }
	    }
	    
	    tempData2 = tempData1[,c(nVar[i], cVar, clustVar)]
	    
	    if (!all(is.na(tempData2[,nVar[i]])) && !all(is.na(tempData2[,cVar]))) {
	      # generates statistics
	      allYStats = NULL
	      yStatEB = NULL
	      yStatEBj = NULL
	      yStatLL = NULL
	      yStatUL = NULL
        yStatHW = NULL
	      if (typeData == "raw") {
	        if (!barsClust || (descStat != "mean")) errBars = FALSE
          for (j in 1:nlevels(tempData2[,clustVar])) {
            #create data subset for clustVar level j
            if (!is.null(clustVar)) {
              tempData = tempData2[which(tempData2[,clustVar] == levels(tempData2[,clustVar])[j]),]
            } else tempData = tempData2
            
	          yStat = NULL
	          if (descStat == "mean") {
	            yStat = tapply(tempData[,nVar[i]], tempData[,cVar], mean, na.rm = TRUE)
	            
	            if (errBars) {
	              ySum <- tapply((tempData[,nVar[i]]), tempData[,cVar], sum, na.rm = TRUE)
	              yFreq <- ySum/yStat
	              yStdDev <- tapply(tempData[,nVar[i]], tempData[,cVar], sd, na.rm = TRUE)
	              yStdErr <- yStdDev / sqrt(yFreq) 
	              
	              if (typeErrBar == "stdDev") {
                  halfWidth = errMult * yStdDev
	              } else if (typeErrBar == "stdErr") {
	                halfWidth = errMult * yStdErr
	              } else if (typeErrBar == "confInt") {
	                halfWidth = qnorm((1+confLev)/2) * yStdErr
	              }
	              lLimit = yStat - halfWidth
	              uLimit = yStat + halfWidth
	              yStatLL = rbind(yStatLL, lLimit)
	              yStatUL = rbind(yStatUL, uLimit)
                yStatHW = rbind(yStatHW, halfWidth)                
	            } # end of if errBars
	          } else if (descStat == "median") {
	            yStat <- tapply(tempData[,nVar[i]], tempData[,cVar], median, na.rm = TRUE)
	          } else if (descStat == "sum") {
	            yStat <- tapply(tempData[,nVar[i]], tempData[,cVar], sum, na.rm = TRUE)
	          } else if (descStat == "freq") {
	            yMean = tapply(tempData[,nVar[i]], tempData[,cVar], mean, na.rm = TRUE)
	            ySum <- tapply((tempData[,nVar[i]]), tempData[,cVar], sum, na.rm = TRUE)
	            yStat <- ySum/yMean
	          }
	          allYStats = rbind(allYStats, yStat)
	        } # end of for (j in 1:nlevels(tempData2[,clustVar]))
	        allYStats = t(allYStats)
          colnames(allYStats) = NULL
	      } else if (typeData == "sumStat") {
	        tempData = na.omit(tempData2[order(tempData2[,cVar]),])
	        allYStats = as.matrix(tempData[,nVar[i]])
	        errBars = FALSE
	      }
	      
	      if (nlevels(tempData[,clustVar]) > 1) { namesX = levels(tempData[,clustVar])
	      } else namesX = NULL
	              
        #sets limits for the numerical variable
	      if (is.na(minValue[i]) && is.na(maxValue[i]) && (!errBars)) {
	        yMinLim = NULL
	        yMaxLim = NULL
	      } else {
	        yMinLim = if(!is.na(minValue[i])) { minValue[i]
	        } else {
	          if(!is.null(yStatLL)) min(0, allYStats, yStatLL, na.rm = TRUE) 
	          else min(0, min(allYStats, na.rm = TRUE)) 
	        }
	        yMaxLim = if(!is.na(maxValue[i])) { maxValue[i]
	        } else {
	          if (barsClust) {
	            if (!is.null(yStatUL)) max(allYStats, yStatUL, na.rm = TRUE) * 1.10
	            else max(allYStats, na.rm = TRUE) * 1.10
	          } else {
	            maxVals = NULL
	            for (s in 1:ncol(allYStats)) maxVals[s] = sum(allYStats[,s], na.rm = TRUE)
	            max(maxVals, na.rm = TRUE) * 1.10
	          }
	        }        
	      }
	      
	      legLab = NULL
	      if (showLeg) legLab = levels(tempData[,cVar]) 
	      yAxisLabel = NULL
	      xAxisLabel = NULL
	      
	      if (barsHoriz) {
	        xVarLim = c(yMinLim,yMaxLim)
	        yVarLim = NULL
	        if (!is.null(xAxisLab)) yAxisLabel = xAxisLab
	        if (!is.null(yAxisLab[i])) xAxisLabel = yAxisLab[i]
	      } else {
	        yVarLim = c(yMinLim,yMaxLim)
	        xVarLim = NULL
	        if (!is.null(yAxisLab[i])) yAxisLabel = yAxisLab[i] else yAxisLabel = NULL
	        if (!is.null(xAxisLab)) xAxisLabel = xAxisLab else xAxisLabel = NULL
	      }
	      
	      barplot1 <- barplot(allYStats, ylab = yAxisLabel, space = c(0.2,1), main = mTitle,
                            names.arg = namesX, horiz = barsHoriz, 
	                          beside = barsClust, col = barColor, xlim = xVarLim, ylim = yVarLim,
	                          density = barDensity, angle = barLineAngle, legend.text = legLab,
	                          args.legend = list(x = legPos, title = legTitle, ncol = legCol), 
                            las = 1)
        
        #draws a box around the plot
	      if (boxed) box(lty = 1)
	      
	      #adds labels per level of the categorical variable axis of the plot
	      if (nlevels(tempData[,clustVar]) > 1 && barsClust && showCatVarLevels) {
	        if (barsHoriz) {
	          mtext(side = 2, at = barplot1, text = rep(levels(tempData[,cVar]), nlevels(tempData[,clustVar])),
	                line = 0, cex = .75, las = axisLabelStyle)
	        } else {
	          mtext(side = 1, at = barplot1, text = rep(levels(tempData[,cVar]), nlevels(tempData[,clustVar])),
	                line = 0, cex = .75, las = axisLabelStyle)
	        } 
	      }										
	      
        #adds subtitle, if any
        if (!is.null(subTitle))
          mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
        
	      #adds a label for the x-axis of the plot, as defined by the user
 	      mtext(side = 1, text = xAxisLabel, line = 2, cex = 0.85)

        #adds error bars to the plot - for barsClust = TRUE and descStat = "mean" only
	      if (errBars) {
	        lLimits = matrix(t(yStatLL),c(nlevels(tempData[,cVar]), nlevels(tempData[,clustVar])))
	        uLimits = matrix(t(yStatUL),c(nlevels(tempData[,cVar]), nlevels(tempData[,clustVar])))
	        if (barsHoriz) { arrows(lLimits, barplot1,  uLimits, barplot1, angle = 90, code = 3, length = 0.1)
	        } else arrows(barplot1, lLimits, barplot1, uLimits, angle = 90, code = 3, length = 0.1)
	      }
	    } #end of if (!all(is.na(tempData[,nVar])) && !all(is.na(tempData[,cVar]))) {
	    
	    #increments file number
	    if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
	      dev.off()
	      k = k + 1
	    }
	    graphNum = graphNum + 1
	  } # end of for (i in 1:length(nVar))

	} # end of for (m in 1:numGroups)
	
}