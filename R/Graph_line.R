# --------------------------------------------------------
# GraphLine: creates a line chart of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# yVars - names of the variable(s) to be plotted against the factor variable
# xVar - name of the factor variable
# lineVars - name of grouping variable(s) for creating different lines
# onePlot - logical; for multiple yVars, whether the yVars will plotted 
#             in the same graph or not
# mTitle - main title for line graph, if any
# yAxisLab - labels to be used for the y-axis variables
# xAxisLab - label to be used for the x-axis
# yMinValue - minimum value(s) for the the y-axis variables
# yMaxValue - maximum value for the y-axis
# axisLabelStyle - style for the axes labels
# byVar - name of categorical variable for which separate graphs are created
# errBars - logical; whether error bars are displayed or not
# typeErrBar - whether error bars displayed (if any) are based on:
#   				standard error or confidence interval
# errMult - multiplier to be used on standard error 
# confLev - confidence level, if confidence interval type error bars are to be displayed
# plotCol - vector of rgb values for the color(s) of the lines (or names of colors)
# showLineLabels - logical; whether the labels beside the lines are displayed or not
# showLeg - logical; whether the graph legend is displayed or not
# legPos - position of the legend (if showLeg is TRUE)
# legTitle - text for the title of the legend, if displayed
# legCol - number of columns used for legend text
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphLine <- function(data, outputPath, yVars, xVar = NULL, lineVars = NULL, mTitle = NULL, yAxisLab = NULL, 
                      xAxisLab = NULL, yMinValue = NULL, yMaxValue = NULL, axisLabelStyle = 1,
                      byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdErr", "confInt"), errMult = 1, confLev = NULL, 
                      plotCol = NULL, showLineLabels = FALSE, showLeg = FALSE, 
                      legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center"), 
                      legTitle = NULL, legCol= 1, boxed = TRUE, linePtTypes, lineTypes,  lineWidths, pointChars, pointCharSizes, 
                      multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, orientGraphs = c("left-right", "top-bottom"))
  UseMethod("GraphLine")

GraphLine.default <- function(data, outputPath, yVars, xVar = NULL, lineVars = NULL, mTitle = NULL, yAxisLab = NULL, 
                      xAxisLab = NULL, yMinValue = NULL, yMaxValue = NULL, axisLabelStyle = 1,
                      byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdErr", "confInt"), errMult = 1, confLev = NULL, 
                      plotCol = NULL, showLineLabels = FALSE, showLeg = FALSE, 
                      legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",	"center"), 
                      legTitle = NULL, legCol= 1, boxed = TRUE, linePtTypes, lineTypes,  lineWidths, pointChars, pointCharSizes, 
                      multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, orientGraphs = c("left-right", "top-bottom")) { 
 
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  if (!any(is.na(as.numeric(levels(data[,xVar]))))) { xNumeric = TRUE
  } else xNumeric = FALSE
  
  #creates a grouping variable (lineVars) if it is not declared
  lineVar = NULL
  if (is.null(lineVars)) { 
    lineVar <- make.unique(c(names(data), "lineVar"))[length(make.unique(c(names(data), "grp.Var")))]
    data[,lineVar] = rep(1,nrow(data))
  } else {
    #combines levels of 2 or more line variables
    if (length(lineVars) == 2) { 
      lineVar = paste(lineVars[1], lineVars[2], sep = "-")
      data[,lineVar] = paste(data[,lineVars[1]], data[,lineVars[2]], sep ="-")
    } else if (length(lineVars) == 3) { 
      lineVar = paste(lineVars[1], lineVars[2], lineVars[3], sep = "-")
      data[,lineVar] = paste(data[,lineVars[1]], data[,lineVars[2]], data[,lineVars[3]], sep ="-")
    } else {
      lineVar = lineVars[1]
      data[,lineVar] = data[,lineVars[1]]
    }  
  }
  
  #converts to factor the grouping variable(s)
  data[,xVar] <- factor(data[,xVar])
  data[,lineVar] <- factor(data[,lineVar]) 
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
    numGraphs = nlevels(data[,byVar])*length(yVars)
  } else numGraphs = length(yVars)
  
  graphNum = 1
  #counts the number of files to save
  k = 1
    
  for (m in 1:numGroups) {
    #creates data by subgroup, if a grouping variable is defined
    if (!is.null(byVar)) {
      tempData1 = data[which(data[,byVar] == levels(data[,byVar])[m]),]
      subTitle = paste(byVar," = ",levels(tempData1[,byVar])[m], sep="")
    } else {
      tempData1 = data
      subTitle = NULL
    }

    for (i in 1:length(yVars)) {
      #creates device for saving graph(s)
      if (!multGraphs) {
        png(filename = paste(outputPath,"linegraph_",k,".png",sep=""))
        par(mfrow=c(numRowsGraphs, numColsGraphs))
        
      } else {
        if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)  {
          widthAdj = numColsGraphs * 480
          heightAdj = numRowsGraphs * 480
          
          png(filename = paste(outputPath,"linegraph_",k,".png", sep=""), width = widthAdj, height = heightAdj)
          
          if (orientGraphs == "top-bottom") {
            par(mfcol=c(numRowsGraphs, numColsGraphs))
          } else if (orientGraphs == "left-right") {
            par(mfrow=c(numRowsGraphs, numColsGraphs))
          }
        }
      }
      
      tempData2 = tempData1[,c(yVars[i], xVar, lineVar)]

      namesX = NULL
        
      if (!all(is.na(tempData2[,yVars[i]])) && !all(is.na(tempData2[,xVar]))) {

        if (xNumeric) { xVals = as.numeric(levels(tempData2[,xVar]))
        } else xVals = c(1:nlevels(tempData2[,xVar]))
          
        # generates statistics
        allYStats = NULL
        yStatEB = NULL
        yStatEBj = NULL
        yStatHW = NULL
        yStatLL = NULL
        yStatUL = NULL 
        yStat = NULL
          
        #determines number of levels of lineVar
        lenlineVar = nlevels(tempData2[,lineVar])
        
        for (j in 1:lenlineVar) {
            
          #create data subset for line j
          if (!is.null(lineVar)) {
            tempData = tempData2[which(tempData2[,lineVar] == levels(tempData2[,lineVar])[j]),]
          } else tempData = tempData2
          
          yStat = tapply(tempData[,yVars[i]], tempData[,xVar], mean, na.rm = TRUE)
              
          if (errBars) {
            ySum = tapply((tempData[,yVars[i]]), tempData[,xVar], sum, na.rm = TRUE)
            yFreq = ySum/yStat
            yStdDev = tapply(tempData[,yVars[i]], tempData[,xVar], sd, na.rm = TRUE)
            yStdErr = yStdDev / sqrt(yFreq) 
            
            if (typeErrBar == "stdErr") {
              halfWidth = errMult * yStdErr
            } else if (typeErrBar == "confInt") {
              halfWidth = qnorm((1+confLev)/2) * yStdErr
            }
            lLimit = yStat - halfWidth
            uLimit = yStat + halfWidth
            yStatLL = rbind(yStatLL, lLimit)
            yStatUL = rbind(yStatUL, uLimit)
            yStatHW = rbind(yStatHW, halfWidth)
          } # end of if (errBars)
          
            allYStats = rbind(allYStats, yStat)
            
        } #end of for (j in 1:nlevels(tempData2[,lineVars]))
        allYStats = t(allYStats)
        allYStats1 = t(allYStats)
        if (!xNumeric) namesX = rownames(allYStats)
        
        #determines lower and upper limits for the y-axis
        yMinLim = if(!is.na(yMinValue[i])) { yMinValue[i]
        } else {
          if(!is.null(yStatLL)) min(allYStats1, yStatLL, na.rm=TRUE) 
          else min(allYStats1, na.rm=TRUE)
        }
        yMaxLim = if(!is.na(yMaxValue[i])) { yMaxValue
        } else {
          if (!is.null(yStatUL)) max(allYStats1, yStatUL, na.rm=TRUE) 
          else max(allYStats1, na.rm=TRUE) 
        }
        
        yVarLim = c(yMinLim, yMaxLim)
        
        legLab = NULL
        q = 1
        
        if (showLeg || showLineLabels) {
          for (r in 1:lenlineVar) {
            if (!all(is.nan(allYStats[,r]))) {
              legLab[q] = levels(tempData[,lineVar])[r]
              q = q + 1
            }
          }
        }
          
        if (xNumeric) {
          xVarLim = c(min(as.numeric(tempData[,xVar])), max(as.numeric(tempData[,xVar])))
        } else  xVarLim = NULL
        
        if (showLineLabels || multGraphs)
          par(mar=c(5, 4, 4, 4) + 0.20)
          
        varToPlot = allYStats[,1] 
          
        plot(xVals, varToPlot, ylab = yAxisLab[i], xlab = xAxisLab, main = mTitle, 
         type = "n", axes=FALSE, col = plotCol, ylim = yVarLim)
          
        if (xNumeric) { axis(1, at=xVals, las = axisLabelStyle) 
        } else { axis(1,at=1:length(namesX), labels = namesX, las = axisLabelStyle)}
        
        axis(2, ylim = yVarLim, las = axisLabelStyle)
        
        #adds subtitle, if any
        if (!is.null(subTitle)) 
          mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
          
        if (showLeg) 
          legend(legPos, legend = legLab, title = legTitle, pch = pointChars, lty = lineTypes, col = plotCol, cex = 0.75, ncol= legCol)
          
        #draws a box around the plot
        if (boxed) box(lty = 1)
          
        #plot lines
        for (l in 1:lenlineVar) {
          varToPlot = allYStats[,l] 
          lines(xVals, varToPlot, type = linePtTypes[l], lwd = lineWidths[l], 
            pch = pointChars[l], cex = pointCharSizes[l], col = plotCol[l], lty = lineTypes[l])
        }
          
        #plot line labels
        if (showLineLabels) {
          for (n in 1:lenlineVar) {
            mtext(side = 4, at = allYStats[nrow(allYStats),n], legLab[n], #yVars[l]
              cex = (pointCharSizes[n]-.25), col = plotCol[n], line = 0.3, las = 1)
          }
        }
          
        #adds error bars to the plot
        if (errBars) {
          lLimits = t(matrix(yStatLL,c(lenlineVar,length(xVals))))
          uLimits = t(matrix(yStatUL,c(lenlineVar,length(xVals))))
          for (h in 1:lenlineVar) {
            arrows(xVals, lLimits[,h], xVals, uLimits[,h], angle = 90, code = 3, length = 0.1, col = plotCol[h])
          } 
        } 
          
      } # end of if (!all(is.na(tempData2[,yVars])) && !all(is.na(tempData2[,xVar]))) 
      if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
        dev.off()
        k = k + 1
      }
      graphNum = graphNum + 1
      
    } # end of for (i in 1:length(yVars))
    
  } # end of for (m in 1:numGroups)
  
}  