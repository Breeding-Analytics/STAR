# --------------------------------------------------------
# GraphBoxplot: creates a boxplot of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# nVar - name of the numeric variable
# cVars - name of grouping variable(s) for dividing data
# mTitle - main title for the boxplot, if any
# nAxisLab - label to be used for the axis of the numeric variable
# cAxisLab - label to be used for axis of the categorical variable
# yMinValue - minimum value for the the y-axis
# yMaxValue - maximum value for the y-axis
# axisLabelStyle - style for the axes labels
# plotHoriz - logical; whether plots are drawn horizontally or not
# byVar - name of categorical variable for which separate graphs are created
# boxColor - vector of rgb values for the color of the box
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphBoxplot <- function(data, outputPath, nVar, cVars = NULL, mTitle = NULL, nAxisLab = NULL, cAxisLab = NULL,
                         yMinValue = NULL, yMaxValue = NULL, axisLabelStyle = 1, plotHoriz = FALSE, byVar = NULL, 
                         boxSize = 0.8, boxColor = NULL, boxFillColor = NULL,
                         medLineType = 1, medLineWidth = 3, medColor = NULL, whiskLineType = 2, whiskLineWidth = 1,
                         whiskColor = NULL, outChar = 1, outCharSize = 1, outColor = NULL,
                         boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                         orientGraphs = c("left-right", "top-bottom"))
  UseMethod("GraphBoxplot")

GraphBoxplot.default <- function(data, outputPath, nVar, cVars = NULL, mTitle = NULL, nAxisLab = NULL, cAxisLab = NULL,
                         yMinValue = NULL, yMaxValue = NULL, axisLabelStyle = 1, plotHoriz = FALSE, byVar = NULL, 
                         boxSize = 0.8, boxColor = NULL, boxFillColor = NULL,
                         medLineType = 1, medLineWidth = 3, medColor = NULL, whiskLineType = 2, whiskLineWidth = 1,
                         whiskColor = NULL, outChar = 1, outCharSize = 1, outColor = NULL,
                         boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                         orientGraphs = c("left-right", "top-bottom")) { 
  
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #creates a grouping variable with ones for all rows if no grouping variable is declared
  if (is.null(cVars)) { 
    cVar <- make.unique(c(names(data), "grp.Var"))[length(make.unique(c(names(data), "grp.Var")))]
    data[,cVar] = rep(1,nrow(data))
  }  
  
  #creates a grouping variable (cVars) if it is not declared
  if (is.null(cVars)) { 
    cVar[1] <- make.unique(c(names(data), "cVars"))[length(make.unique(c(names(data), "grp.Var")))]
    data[,cVar[1]] = rep(1,nrow(data))
  } else {
    #combines levels of 2 or more line variables
    if (length(cVars) == 2) { 
      cVar = paste(cVars[1], cVars[2], sep = "-")
      data[,cVar] = paste(data[,cVars[1]], data[,cVars[2]], sep ="-")
    } else if (length(cVars) == 3) { 
      cVar = paste(cVars[1], cVars[2], cVars[3], sep = "-")
      data[,cVar] = paste(data[,cVars[1]], data[,cVars[2]], data[,cVars[3]], sep ="-")
    } else {
      cVar = cVars[1]
      data[,cVar] = data[,cVars[1]]
    }
  }
  
  #converts to factor the grouping variable(s)
  data[,cVar] <- factor(data[,cVar]) 
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
    
    tempData1 = tempData1[,c(nVar, cVar)]
    
    if (all(is.na(tempData1[,cVar]))) {
      tempData1[,cVar] = rep(1,nrow(tempData1))
    }  
    
    tempData1 = tempData1[order(tempData1[,cVar]),]
    
    namesX = levels(tempData1[,cVar])

    for (i in 1:length(nVar)) {
      #creates device for saving graph(s)
      if (!multGraphs) {
        png(filename = paste(outputPath,"boxplot",k,".png",sep=""))
        par(mfrow=c(numRowsGraphs, numColsGraphs))
        
      } else {
        if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)  {
          widthAdj = numColsGraphs * 480
          heightAdj = numRowsGraphs * 480
          png(filename = paste(outputPath,"boxplot",k,".png",sep=""), width = widthAdj, height = heightAdj)
          
          if (orientGraphs == "top-bottom") {
            par(mfcol=c(numRowsGraphs, numColsGraphs))
          } else if (orientGraphs == "left-right") {
            par(mfrow=c(numRowsGraphs, numColsGraphs))
          }
        }
      }
      
      tempData = na.omit(tempData1[,c(nVar[i],cVar)])
      if (nrow(tempData) > 0) {
        
        #determines lower and upper limits for the y-axis
        if (is.na(yMinValue[i]) && is.na(yMaxValue[i])) {
          yVarLim = NULL
        }
        else {
          yMinLim = if(!is.na(yMinValue[i])) { yMinValue[i]
          } else min(tempData[,nVar[i]], na.rm = TRUE) 
          
          yMaxLim = if(!is.na(yMaxValue[i])) { yMaxValue[i]
          } else max(tempData[,nVar[i]], na.rm = TRUE) 
          
          yVarLim = c(yMinLim,yMaxLim)
        }
        
        if (plotHoriz) {
          yAxisLab = cAxisLab
          xAxisLab = nAxisLab[i]
        } else {
          xAxisLab = cAxisLab
          yAxisLab = nAxisLab[i]
        }
        
        boxplot(tempData[,nVar[i]] ~ tempData[,cVar], ylab = yAxisLab, xlab = xAxisLab, main = mTitle, 
                ylim = yVarLim, las = axisLabelStyle, frame.plot = boxed, horizontal = plotHoriz, names.arg = namesX, 
                pars = list(boxwex = boxSize, boxcol = boxColor, 
                boxfill = boxFillColor, medlty = medLineType, medlwd = medLineWidth, medcol = medColor,
                whisklty = whiskLineType, whisklwd = whiskLineWidth, whiskcol = whiskColor, staplecol = whiskColor,
                outpch = outChar, outcex = outCharSize, outcol = outColor))
        
        #adds subtitle, if any
        if (!is.null(subTitle))
          mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
        
      }
      
      #increments file number
      if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
        dev.off()
        k = k + 1
      }
      graphNum = graphNum + 1
      
    }  #end of if (!all(is.na(tempData[,nVar])) && !all(is.na(tempData[,cVar]))) {    

  }
  
}