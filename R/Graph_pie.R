# --------------------------------------------------------
# GraphPie: creates a piechart of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# cVar - name of the categorical variable
# mTitle - main title for bar graph, if any
# wVar - name of weight (frequencies) variable, if any
# byVar - name of categorical variable for which separate graphs are created
# sliceColors - colors for the slices
# dispLabels - whether labels are used to describe the slices
# useLabels - type of label to use
# showLeg - logical; whether the graph legend is displayed or not
# legPos - position of the legend (if showLeg is TRUE)
# legTitle - text for the title of the legend, if displayed
# legCol - number of columns used for legend text
# pieDensity - line density per slice 
# pieLineAngle - angles for the lines per slice
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphPie <- function(data, outputPath, cVar, mTitle = NULL, wVar = NULL,
                     byVar = NULL, sliceColors, dispLabels = FALSE, useLabels = c("counts", "pcts", "cats", "pctcats", "countcats"), 
                     showLeg = FALSE,legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright",
                                                "right", "center"), legTitle = NULL, legCol = 1, pieDensity, pieLineAngle, boxed = TRUE,
                     multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1,
                     orientGraphs = c("left-right", "top-bottom"))
  UseMethod("GraphPie")

GraphPie.default <- function(data, outputPath, cVar, mTitle = NULL, wVar = NULL,
								byVar = NULL, sliceColors, dispLabels = FALSE, useLabels = c("counts", "pcts", "cats", "pctcats", "countcats"), 
                showLeg = FALSE,legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright",
                    "right", "center"), legTitle = NULL, legCol = 1, pieDensity, pieLineAngle, boxed = TRUE,
								multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1,
								orientGraphs = c("left-right", "top-bottom")) { 
	
	#reads data
	if (is.character(data)) { data <- eval(parse(text = data)) }
	
	#converts to factor the grouping variable(s)
	if (!is.null(cVar)) { data[,cVar] <- factor(data[,cVar]) }
	if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }	
	
	if (!multGraphs) {
		numRowsGraphs = 1
		numColsGraphs = 1
	} 
	
	#determines number of cells allocated for graphs (esp. for multiple graphs)
 	numCells = numRowsGraphs * numColsGraphs
	
	#determines number of graphs to be created
	if (!is.null(byVar)) {
		numGraphs = (nlevels(data[,byVar]))
	} else numGraphs = 1

 	graphNum = 1

  #counts the number of files to save
	k = 1

	for (m in 1:numGraphs) {
		#creates device for saving graph(s)
		if (!multGraphs) {
			png(filename = paste(outputPath,"piechart",k,".png",sep=""))
			par(mfrow=c(numRowsGraphs, numColsGraphs))
			
		} else {
				if (graphNum == 1 || graphNum %% numCells == 1 || numCells == 1)	{
				  widthAdj = numColsGraphs * 480
				  heightAdj = numRowsGraphs * 480
				  
					png(filename = paste(outputPath,"piechart",k,".png", sep=""), width = widthAdj, height = heightAdj)
				
					if (orientGraphs == "top-bottom") {
						par(mfcol=c(numRowsGraphs, numColsGraphs))
					} else if (orientGraphs == "left-right") {
						par(mfrow=c(numRowsGraphs, numColsGraphs))
					}
				}
		}
		
		#creates data by subgroup, if a grouping variable is defined
		if (!is.null(byVar)) {
			tempData = data[which(data[,byVar] == levels(data[,byVar])[m]),]
			subTitle = paste(byVar,"=",levels(tempData[,byVar])[m], sep=" ")
		} else {
			tempData = data
			subTitle = NULL
		}

		if (length(na.omit(tempData[,cVar])) > 0) {
		  # generates counts
		  
		  if (is.null(wVar)) {
		    countsData = table(tempData[,cVar])
		  } else {
		    counts1 = tempData[order(tempData[,cVar]),]
		    countsData = as.table(counts1[,wVar]) 
		    names(countsData) = levels(data[,cVar])
		  }
		  
		  props1 = prop.table(countsData)
		  propData = as.data.frame(props1)
		  propData$Freq = propData$Freq * 100
		  colnames(propData) = c("Var1", "Pcts")
		  pieData = merge(countsData, propData, by = "Var1" )
		  pieData$Pcts = format(pieData$Pcts, digits = 4)
		  namesX = levels(tempData[,cVar])
		  
		  legLab = NULL
		  if (showLeg) legLab = namesX
		  
		  labelText = NULL
		  
		  if (dispLabels) {
		    if (useLabels == "counts") labelText = pieData$Freq
		    else if (useLabels == "pcts") labelText = paste(pieData$Pcts,"%", sep="")  
		    else if (useLabels == "cats") labelText = pieData$Var1  
		    else if (useLabels == "pctcats") labelText = paste(pieData$Var1, " ", pieData$Pcts, "%", sep = "")
		    else if (useLabels == "countcats") labelText = paste(pieData$Var1, pieData$Freq)
		  } else labelText = NA
		  
		  pie(pieData$Freq, main = mTitle, labels = labelText, col = sliceColors,
		      density = pieDensity, angle = pieLineAngle, cex = 0.8, border = "black")
		  	  
		  #adds subtitle, if any
		  if (!is.null(subTitle)) 
		    mtext(side = 3, text = subTitle, line = 0.25, cex = 0.9)
		  
		  #draws a box around the plot
		  if (boxed) box(lty = 1)
		  
		  if (showLeg) 
		    legend(legPos, legend = legLab, title = legTitle, ncol = legCol, fill = sliceColors, density = pieDensity * 2, angle = pieLineAngle, cex = 0.8)
		  
		  #increments file number
		  if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
		    dev.off()
		    k = k + 1
		  }
		  
		  graphNum = graphNum + 1
		}
    
	}
	
}