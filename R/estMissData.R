# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# estMissData: Functions for Estimating Missing Data from ANOVA
# Created by: Alaine A. Gulles 07.05.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.27.2012
# -------------------------------------------------------------------------------

estMissData <- function(design, data, respvar, factor1, factor2, factor3, factor4, rep1, rep2) UseMethod("estMissData")

estMissData.default <- function(design, data, respvar, factor1, factor2, factor3, factor4, rep1, rep2) {
	availableDesign <- c("CRD", "RCBD", "LSD", "SplitCRD", "SplitRCBD", "SplitLSD", "Strip", "Split2CRD", "Split2RCBD", "Split2LSD", "Strip-Split",	"Split3CRD", "Split3RCBD", "Split3LSD", "Strip-Split2")
	design <- availableDesign[match(design, availableDesign)]
	designChoice <- match(design, availableDesign)
    	if (is.character(data)) { data <- eval(parse(text = data)) }
	tempData <- data
	blk <- c(rep1, rep2)
	factor <- c(factor1, factor2, factor3, factor4)

     	switch(designChoice,
            {modelRHS <- paste(paste(factor1, collapse = "*", sep = ""))},
            {modelRHS <- paste(rep1, " + ", paste(factor1, collapse = "*", sep = ""), sep = "")},
            {modelRHS <- paste(rep1, " + ", rep2, " + ", paste(factor1, collapse = "*", sep = ""), sep = "")},
            {modelRHS <- paste(paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = "") , ")/",  rep1, ")", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error(", paste(c(rep1, factor1), collapse = ":", sep = ""), "/(", paste(factor1, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error(", paste(c(rep1,rep2, factor1), collapse = ":", sep = ""), "/(", paste(factor1, collapse = "*", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2), collapse = "*", sep = ""), " + Error((", paste(c(rep1,factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", rep2 ," + ",paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(c(factor1, factor2), collapse = ":", sep = ""),"))", sep = "")},
            {modelRHS <- paste(paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", rep1, ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", rep1, "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""),         "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""),         "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(factor3, collapse = ":", sep = ""), "))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", rep2 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, rep2, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, rep2, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(factor3, collapse = ":", sep = ""),"))", sep = "")},
            {modelRHS <- paste(rep1 ," + ", paste(c(factor1, factor2, factor3, factor4), collapse = "*", sep = ""), " + Error((", paste(c(rep1, factor1), collapse = ":", sep = ""), "/", paste(factor1, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor2), collapse = ":", sep = ""), "/", paste(factor2, collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2), collapse = ":", sep = ""), "/", paste(c(factor1, factor2), collapse = ":", sep = ""), ") + (", paste(c(rep1, factor1, factor2, factor3), collapse = ":", sep = ""), "/", paste(c(factor3), collapse = ":", sep = ""), "))", sep = "")}
     	)	


	if (design == "LSD" || design == "SplitLSD" || design == "Split2LSD" || design == "Split3LSD") {
		# obtain Row means
		RowMeans <- data.frame(levels(tempData[,blk[1]]), tapply(tempData[,respvar], tempData[,blk[1]], mean, na.rm = TRUE))
		colnames(RowMeans) <- c(blk[1], paste(blk[1], "means", sep = ""))
		# obtain Column means
		ColMeans <- data.frame(levels(tempData[,blk[2]]), tapply(tempData[,respvar], tempData[,blk[2]], mean, na.rm = TRUE))
		colnames(ColMeans) <- c(blk[2], paste(blk[2], "means", sep = ""))
		# merge the data
		newData <- data.frame(merge(merge(tempData, RowMeans, by = blk[1]), ColMeans, by = blk[2]), remarks = "obs")
		newData$remarks <- as.character(newData$remarks)
		newData[is.na(newData[,respvar]),"remarks"] <- "estimate"
		newData$Means <- rowMeans(newData[,match(c(paste(blk[1], "means", sep = ""),paste(blk[2], "means", sep = "")), names(newData))])
		newData[is.na(newData[,respvar]),respvar] <- newData[is.na(newData[,respvar]),paste("Means", sep = "")]
		newData <- newData[,-I(match(c(paste(blk[1], "means", sep = ""), paste(blk[2], "means", sep = "")),names(newData)))]
		colnames(newData)[match("Means", names(newData))] <- "tempSum"
	} else {
		# obtain the means per replicate/blk
		RepMeans <- data.frame(levels(tempData[,blk]), tapply(tempData[,respvar], tempData[,blk], mean, na.rm = TRUE))
		colnames(RepMeans) <- c(blk, paste(respvar, "means", sep = ""))
		newData <- data.frame(merge(tempData, RepMeans, by = blk), remarks = "obs")
		newData$remarks <- as.character(newData$remarks)
		newData[is.na(newData[,respvar]),"remarks"] <- "estimate"
		newData[is.na(newData[,respvar]),respvar] <- newData[is.na(newData[,respvar]),paste(respvar, "means", sep = "")]
		colnames(newData)[match(paste(respvar, "means", sep = ""),names(newData))] <- "tempSum"
		#newData <- newData[,-I(match(paste(respvar, "means", sep = ""),names(newData)))]
		#if (design == "CRD" || design == "SplitCRD" || design == "Split2CRD" || design == "Split3CRD") { 
		#	myformula <- paste(respvar, " ~ ", paste(factor, collapse = ":", sep = ""), sep = "")
		#} else { myformula <- paste(respvar, " ~ ", paste(blk, collapse = " + ", sep = "")," + ", paste(factor, collapse = ":", sep = ""), sep = "") }
	}
	myformula <- paste(respvar, " ~ ", modelRHS, sep = "")
	stable <- FALSE 
	iterationCounter <- 0
	while(!stable || iterationCounter > 100) {
    		iterationCounter <- iterationCounter + 1
		result <- suppressWarnings(aov(formula(myformula), data = newData))
		newData$predval <- PredictedValues(result)
		estimatedData <- subset(newData, remarks == "estimate")
		if (all(abs(estimatedData[,respvar]-estimatedData[,"predval"])/estimatedData[,respvar] < 0.001)) { 
			stable <- TRUE
		} 
    		newData[newData[,"remarks"] == "estimate",respvar] <- newData[newData[,"remarks"] == "estimate","predval"]
    		newData[,"tempSum"] <- newData[,"tempSum"] + newData[,"predval"]
    		newData <- newData[,-I(match("predval", names(newData)))]
	}
	if (!stable) { 
    		newData[,"tempSum"] <- newData[,"tempSum"]/iterationCounter
    		newData[newData[,"remarks"] == "estimate",respvar] <- newData[newData[,"remarks"] == "estimate","tempSum"]
  	}
	newData <- newData[,-I(match("tempSum", names(newData)))]
	newVarName <- make.unique(c(colnames(newData), paste(respvar,"remarks", sep = "_")), sep = "")
	colnames(newData)[ncol(newData)] <- newVarName[length(newVarName)]
	return(newData)
} ## end function

