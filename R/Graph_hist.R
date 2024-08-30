# --------------------------------------------------------
# GraphHist: creates a histogram of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# nVar - name of the numeric variable
# mTitle - main title for the boxplot, if any
# yAxisLabs - vector of labels to be used for the y-axis
# xAxisLabs - vector of labels to be used for the x-axis 
# yMinValues - vector of minimum values for the the y-axis
# yMaxValues - vector of maximum values for the y-axis
# xMinValues - vector of minimum values for the the x-axis
# xMaxValues - vector of maximum values for the x-axis
# useFreq - logical; whether to use frequencies in the y-axis or not (probabilities)
# numBins - number of bins; by default, it is determined using Sturges' formula
# axisLabelStyle - style for the axes labels
# byVar - name of categorical variable for which separate graphs are created
# barColor - rgb values for the color of the bars
# dispCurve - logical; whether the kernel density estimates are displayed (in a curve) or not
# lineType - number for types of lines used
# lineColor - rgb values for the color of the line, if any
# lineWidth - width for the line, if any
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphHist <- function(data, outputPath, nVar, mTitle = NULL, yAxisLabs, xAxisLabs, 
                      yMinValues, yMaxValues, xMinValues, xMaxValues, 
                      useFreq = TRUE, numBins = "Sturges", axisLabelStyle = 1, byVar = NULL, 
                      barColor, dispCurve = FALSE, lineType, lineColor, lineWidth,
                      boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                      orientGraphs = c("left-right", "top-bottom")) 
  UseMethod("GraphHist")

GraphHist.default <- function(data, outputPath, nVar, mTitle = NULL, yAxisLabs, xAxisLabs, 
              yMinValues, yMaxValues, xMinValues, xMaxValues, 
              useFreq = TRUE, numBins = "Sturges", axisLabelStyle = 1, byVar = NULL, 
              barColor, dispCurve = FALSE, lineType, lineColor, lineWidth,
              boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
              orientGraphs = c("left-right", "top-bottom")) { 
  
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #converts to factor the grouping variable(s)
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
    numGraphs = nlevels(data[,byVar]) * length(nVar)
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

    for (j in 1:length(nVar)) {
      if (all(is.na(tempData1[,nVar[j]]))) {
        nVar[j] <- NA
      }
    }
    nVar = na.omit(nVar)
    
    if (length(nVar == 1)) { 
      cVar <- make.unique(c(names(tempData1), "grp.Var"))[length(make.unique(c(names(tempData1), "grp.Var")))]
      tempData1[,cVar] = c(1:nrow(tempData1))
    }
    tempData = tempData1[,c(nVar,cVar)]
    
    if (!all(is.na(tempData[,nVar]))) { 
      
      histData = NULL  
      minCount = NULL
      maxCount = NULL
      minDensity = NULL
      maxDensity = NULL
      
      for (i in 1:length(nVar)) {

        #creates device for saving graph(s)
        if (!multGraphs) {
          png(filename = paste(outputPath,"histogram", k,".png",sep=""))
          par(mfrow=c(numRowsGraphs, numColsGraphs))
                
        } else {
          if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)  {
            widthAdj = numColsGraphs * 480
            heightAdj = numRowsGraphs * 480
            
            png(filename = paste(outputPath,"histogram", k,".png", sep=""), width = widthAdj, height = heightAdj)
             
            if (orientGraphs == "top-bottom") {
              par(mfcol=c(numRowsGraphs, numColsGraphs))
            } else if (orientGraphs == "left-right") {
              par(mfrow=c(numRowsGraphs, numColsGraphs))
            }
          }
        }
        
        histData = hist(tempData[,nVar[i]], breaks = numBins, plot = F)  
        minCount = min(histData$counts, na.rm = TRUE)
        maxCount = max(histData$counts, na.rm = TRUE)
        minDensity = min(histData$density, na.rm = TRUE)
        maxDensity = max(histData$density, na.rm = TRUE)

        #determines lower and upper limits for the x-axis
        xMinLim = if(!is.na(xMinValues[i])) { xMinValues[i]
        } else min(0, min(tempData[,nVar[i]], na.rm = TRUE)) 
        xMaxLim = if(!is.na(xMaxValues[i])) { xMaxValues[i]
        } else max(tempData[,nVar[i]], na.rm = TRUE) * 1.25
        
        xVarLim = c(xMinLim,xMaxLim)  
        
        #determines lower and upper limits for the y-axis
        if (is.na(yMinValues[i]) && is.na(yMaxValues[i])) {
          yVarLim = NULL
        } else {
          if (useFreq) {
            yMinLim = if(!is.na(yMinValues[i])) { yMinValues[i]
            } else minCount
            yMaxLim = if(!is.na(yMaxValues[i])) { yMaxValues[i]
            } else {
              maxCount
            }
          } else {
            yMinLim = if(!is.na(yMinValues[i])) { yMinValues[i]
            } else minDensity
            yMaxLim = if(!is.na(yMaxValues[i])) { yMaxValues[i]
            } else {
              if (dispCurve) maxDensity * 1.25
              else maxDensity
            }
          }
          yVarLim = c(yMinLim,yMaxLim)  
        }
      
        plot(histData, xlim = xVarLim, ylim = yVarLim, col = barColor, freq = useFreq, main = mTitle,
         ylab = yAxisLabs[i], xlab = xAxisLabs[i], las = axisLabelStyle)
      
        #adds subtitle, if any
        if (!is.null(subTitle)) 
          mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
        
        if (dispCurve && !useFreq)
          lines(density(tempData[,nVar[i]], na.rm = TRUE), col = lineColor, lty = lineType, lwd = lineWidth)
      
        #draws a box around the plot
        if (boxed) box(lty = 1)
        
        #increments file number
        if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
          dev.off()
          k = k + 1
        }
        graphNum = graphNum + 1
        
      } #end of for (i in 1:length(nVar))

    } # end of if (!all(is.na(tempData[,nVar])))
    
  } # end of for (m in 1:numGroups)
	
}