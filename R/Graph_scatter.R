# --------------------------------------------------------       
# GraphScatter: creates a line chart of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# xVar - list of variables for the x-axis
# yVar - list of variables for the y-axis
# mTitle - main title for bar graph, if any
# xAxisLab - labels to be used for x-axis (for the variables in the x-axis)
# yAxisLab - labels to be used for the y-axis (for the variables in the y-axis)
# xMinValue - minimum value for the the x-axis
# xMaxValue - maximum value for the x-axis
# yMinValue - minimum value for the the y-axis
# yMaxValue - maximum value for the y-axis
# axisLabelStyle - style for the axes labels
# byVar - name of categorical variable for which separate graphs are created
# boxed - logical; whether a box is drawn around the plot or not
# pointCol - rgb values for the color of the points
# pointChar - code (number) for the character used for points
# pointCharSize - size of the character used for points 
# dispLine - logical; whether to display the regression line or not
# lineType - type of line used
# lineCol - color of line
# lineWidth - width of the line
# dispR2P - logical; whether to display R^2 and p-value or not
# r2PPos - position of the legend (if showLeg is TRUE)
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphScatter <- function(data, outputPath, xVar, yVar, mTitle = NULL, xAxisLab = NULL, 
                         yAxisLab = NULL, xMinValue = NULL, xMaxValue = NULL , yMinValue = NULL, yMaxValue = NULL,
                         axisLabelStyle = 1, byVar = NULL, boxed = TRUE, pointCol, pointChar, pointCharSize, 
                         dispLine = FALSE, lineType = NULL, lineCol = NULL, lineWidth, dispR2P = FALSE, 
                         r2PPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center"),
                         multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                         orientGraphs = c("left-right", "top-bottom"))
  UseMethod("GraphScatter")

GraphScatter.default <- function(data, outputPath, xVar, yVar, mTitle = NULL, xAxisLab = NULL, 
		            yAxisLab = NULL, xMinValue = NULL, xMaxValue = NULL , yMinValue = NULL, yMaxValue = NULL,
                axisLabelStyle = 1, byVar = NULL, boxed = TRUE, pointCol, pointChar, pointCharSize, 
                dispLine = FALSE, lineType = NULL, lineCol = NULL, lineWidth, dispR2P = FALSE, 
                r2PPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center"),
                multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                orientGraphs = c("left-right", "top-bottom")) { 
	
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #converts to factor the grouping variable, if any
  if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }  
  
  if (!multGraphs) {
    numRowsGraphs = 1
    numColsGraphs = 1
  } 
  
  #determines number of cells allocated for graphs (esp. for multiple graphs)
  numCells = numRowsGraphs * numColsGraphs
  
  #determines number of graphs to be created
  if (!is.null(byVar)) {
    byNLev = nlevels(data[,byVar])
    numGraphs = byNLev * length(xVar) * length(yVar)
  } else {
    byNLev = 1
    numGraphs = length(xVar) * length(yVar)
  }  
  
  graphNum = 1
  #counts the number of files to save
  k = 1
  
  for (m in 1:byNLev) {

    #creates data by subgroup, if a grouping variable is defined
    if (!is.null(byVar)) {
      tempData1 = data[which(data[,byVar] == levels(data[,byVar])[m]),]
      subTitle = paste(byVar,"=",levels(tempData1[,byVar])[m], sep=" ")
    } else {
      tempData1 = data
      subTitle = NULL
    }

    yStat = NULL
    for (i in 1:length(xVar)) {
      for (j in 1:length(yVar)) {
        #creates device for saving graph(s)
        if (!multGraphs) {
          png(filename = paste(outputPath,"scatter_",k,".png",sep=""))
          par(mfrow=c(numRowsGraphs, numColsGraphs))
          
        } else {
          if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)  {
            widthAdj = numColsGraphs * 480
            heightAdj = numRowsGraphs * 480
            
            png(filename = paste(outputPath,"scatter_",k,".png", sep=""), width = widthAdj, height = heightAdj)
            
            if (orientGraphs == "top-bottom") {
              par(mfcol=c(numRowsGraphs, numColsGraphs))
            } else if (orientGraphs == "left-right") {
              par(mfrow=c(numRowsGraphs, numColsGraphs))
            }
          }
        }
        
        tempData = na.omit(tempData1[,c(xVar[i], yVar[j])])
        if (nrow(tempData) > 0) {
          yMinLim = NULL
          yMaxLim = NULL
          yVarLim = NULL
          xMinLim = NULL
          xMaxLim = NULL
          xVarLim = NULL
           
          #determines lower and upper limits for the y-axis
          if (is.na(yMinValue[j]) && is.na(yMaxValue[j])) {
            yVarLim = NULL
          }
          else {
            yMinLim[j] = if(!is.na(yMinValue[j])) { yMinValue[j]
            } else min(0, min(tempData[,yVar[j]], na.rm = TRUE)) 
            
            yMaxLim[j] = if(!is.na(yMaxValue[j])) { yMaxValue[j]
            } else max(tempData[,yVar[j]], na.rm = TRUE) 

            yVarLim = c(yMinLim[j], yMaxLim[j])
          }
          
          #determines lower and upper limits for the x-axis
          if (is.na(xMinValue[i]) && is.na(xMaxValue[i])) {
            xVarLim = NULL
          }
          else {
            xMinLim[i] = if(!is.na(xMinValue[i])) { xMinValue[i]
            } else min(0, min(tempData[,xVar[i]], na.rm = TRUE)) 
            
            xMaxLim[i] = if(!is.na(xMaxValue[i])) { xMaxValue[i]
            } else max(tempData[,xVar[i]], na.rm = TRUE) 
            
            xVarLim = c(xMinLim[i], xMaxLim[i])
          }
          
          plot(tempData[,xVar[i]], tempData[,yVar[j]], ylab = yAxisLab[j], xlab = xAxisLab[i], main = mTitle, 
               type = "p", pch = pointChar, cex = pointCharSize, col = pointCol, 
               xlim = xVarLim, ylim = yVarLim, las = axisLabelStyle)
          
          #adds subtitle, if any
          if (!is.null(subTitle))
            mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
          
          #draws a box around the plot
          if (boxed) box(lty = 1)
          
          #performs regression
          if (nrow(tempData) > 1) {
            regMod = lm(tempData[,yVar[j]] ~ tempData[,xVar[i]])
            regModSum = summary(regMod)
            
            #overlay regression line
            if (dispLine) {
              abline(regMod, cex = pointCharSize, col = lineCol, lty = lineType)
              
              #prints R2 and p-value
              if (dispR2P) {
                rSqr = regModSum$r.squared * 100
                pValue = pf(regModSum$fstatistic[1],regModSum$fstatistic[2],regModSum$fstatistic[3],lower.tail=FALSE)
                r2P = vector('expression',2)
                r2P[1] = substitute(expression(italic(R)^2~"(%)" == rSqr), list(rSqr = format(rSqr,digits = 3)))[2]
                r2P[2] = substitute(expression(italic(p) == pValue), list(pValue = format(pValue, digits = 4)))[2]
                legend(r2PPos,  r2P, bty = "n", cex = 0.75)          
              }
            }
          }
        }  #end of if(nrow(data) > 0)  
 
        
        #increments file number
        if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
          dev.off()
          k = k + 1
        }
        graphNum = graphNum + 1
      }
    }
  }
}
