# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# GenerateBalanceData: 
# Created by: Alaine A. Gulles 07.18.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.19.2012
# -------------------------------------------------------------------------------

GenerateBalanceData <- function(data, respvar, factor, blk, design) UseMethod("GenerateBalanceData")

GenerateBalanceData.default <- function(data, respvar, factor, blk, design) {
	if (is.character(data)) { data <- eval(parse(text = data)) } 
	if (!is.data.frame(data)) { stop("The object must be a data frame.") }

	availableDesign <- c("CRD", "RCBD", "LSD", "SplitCRD", "SplitRCBD", "SplitLSD", "Strip", "Split2CRD", "Split2RCBD", "Split2LSD", "Strip-Split",	"Split3CRD", "Split3RCBD", "Split3LSD", "Strip-Split2")
	design <- availableDesign[match(design, availableDesign)]
	
	rawData <- data[,sort(match(c(respvar, factor, blk), names(data)))]
	numObs <- replications(paste(respvar, " ~ ", paste(c(factor, blk[1]), collapse = ":"), sep = ""), rawData)
	if (is.list(numObs)) {
		if (max(numObs[[1]]) != 1) { stop(paste("No unique balance data can be generated for variable ", respvar,".", sep = "")) }
		base <- function(myrep, mydata) {
			temp <- replications(paste(respvar, " ~ ", paste(c(factor, myrep), collapse = ":"), sep = ""), rawData)
			fullData <- as.data.frame.table(temp[[1]])
			fullData <- fullData[,-I(ncol(fullData))]
			fullData <- merge(fullData, mydata, by = c(factor, myrep), all.x = TRUE, order = TRUE)
			return(fullData)
		}
		
		if (design == "LSD" || design == "SplitLSD" || design == "Split2LSD" || design == "Split3LSD") {
			rowData <- base(blk[1], mydata = rawData)
			colData <- base(blk[2], mydata = rawData)
			newData <- merge(rowData, colData, by = c(factor, blk, respvar))
			newData <- merge(newData, data, by = c(factor, blk, respvar))
			return(newData[,c(match(names(data), names(newData)))])
		} else { 
			newData <- base(blk[1], mydata = data)
			return(newData[,c(match(names(data), names(newData)))])
		}
	} else { return(data)	}
} ### end stmt -- GenerateBalanceData
