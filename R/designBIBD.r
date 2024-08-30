# -------------------------------------------------------------------------------------
# designBIBD: Generate randomization for Balanced Incomplete Block.
# Created by: Alaine A. Gulles 04.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 09.10.2012
# Note: This is a modified version of the function design.bib in package agricolae 
# -------------------------------------------------------------------------------------

designBIBD <- function(generate, blkSize, trial = 1, numFieldRow = 1, rowPerBlk = 1, serpentine = FALSE, display = TRUE, file = NULL) UseMethod("designBIBD")

designBIBD.default <- function(generate, blkSize, trial = 1, numFieldRow = 1, rowPerBlk = 1, serpentine = FALSE, display = TRUE, file = NULL) {

	if (missing(generate)) { stop("The argument 'generate' is missing.") }
	if (missing(blkSize)) { stop("The argument 'blkSize' is missing.") }
	if (length(generate) > 1) { stop("The argument 'generate' must be of lenght 1.") }

	generate <- FactorList(generate)
	if (blkSize >= length(generate[[1]])) { stop("The argument 'blkSize' should be less than levels of the factor.") }

	tempComb <- GenerateFactor(generate, times = 1)
	trmtComb <- combn(generate[[1]], blkSize)
	numBlks <- ncol(trmtComb)
     
	r <- factorial(length(generate[[1]]) - 1)/(factorial(blkSize - 1)*factorial(length(generate[[1]]) - blkSize))
	lambda <- factorial(length(generate[[1]]) - 2)/(factorial(blkSize - 2)*factorial(length(generate[[1]]) - blkSize))
	randomize <- NULL
     
     plan <- list()
     plotNum <- NULL
     
     if (rowPerBlk == 1) serpentine <- FALSE
     
     numBlkRow <- numFieldRow/rowPerBlk
     numBlkCol <- numBlks/numBlkRow
     colPerBlk <- blkSize/rowPerBlk

	if ((blkSize*numBlks)%%numFieldRow != 0) { stop("Total number of plots should be divisible by the number of field rows.") }
	if ((blkSize%%rowPerBlk) != 0) { stop("Total number of plots per block should be divisible by the number of rows within block") }
	
     
	for (i in (1:trial)) {
		temp <- trmtComb[,sample(1:numBlks, numBlks)]
          plan[[i]] <- matrix(0, nrow = numFieldRow, ncol = numBlkCol*colPerBlk)
          if (i == 1) plotNum <- plan[[i]]
          for (j in (1:numBlks)) { 
               temp[,j] <- temp[sample(1:blkSize,blkSize),j]
               tempPlan <- matrix(temp[,j], nrow = rowPerBlk,ncol = blkSize/rowPerBlk, byrow =TRUE)
               plotLabel <- as.numeric(paste(j, paste(c(rep(0, max(nchar(nrow(temp))))), collapse = ""), sep = ""))+1:nrow(temp)
               tempPlotNum <- matrix(plotLabel, nrow = rowPerBlk,ncol = blkSize/rowPerBlk, byrow =TRUE)
               if (serpentine) { for (k in seq(2, rowPerBlk, by = 2)) { tempPlotNum[k, ] <- rev(tempPlotNum[k, ]) }}
               if (j%%numBlkCol != 0) { colIndex <- j%%numBlkCol } else { colIndex <- numBlkCol }
               if (colIndex == 1) { colIndexLL <- colIndex  } else { colIndexLL <- colIndexLL + colPerBlk }
               colIndexUL <- colIndexLL + colPerBlk - 1
               rowIndex <- ceiling(j/numBlkCol)
               rowIndexLL <- (rowIndex * rowPerBlk) - rowPerBlk + 1
               rowIndexUL <- rowIndexLL + rowPerBlk - 1 
               plan[[i]][rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlan
               if (i == 1) plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlotNum
               tempFieldOrder <- merge(as.data.frame.table(matrix(rowIndexLL:rowIndexUL,nrow = rowPerBlk, ncol = colPerBlk, byrow = FALSE)), 
                                       as.data.frame.table(matrix(colIndexLL:colIndexUL,nrow = rowPerBlk, ncol = colPerBlk, byrow = TRUE)), 
                                       by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
               tempOrder <- merge(as.data.frame.table(tempPlan), as.data.frame.table(tempPlotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
               randomize <- rbind(randomize, merge(data.frame(Trial = i, Blocks = j, Trmt = temp[,j]),
                                                   merge(tempOrder, tempFieldOrder, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))[,3:6], by.x = "Trmt", by.y = "Freq.x.x"))
               
		}
		dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
		if (i == 1) dimnames(plotNum) <- dimnames(plan[[i]])
	}
	names(plan) <- paste("Trial", 1:trial, sep = "")
	names(randomize)[(ncol(randomize)-2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
	randomize <- randomize[, c("Trial", "Blocks", "Trmt", "PlotNum", "FieldRow", "FieldCol")]
	names(randomize)[3] <- names(generate)[1]
	randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$PlotNum),] # fieldbook is sorted by Trial, then by Block then by PlotNum
	rownames(randomize) <- 1:nrow(randomize)
     

	if (display) {
		cat(toupper("Design Properties:"),"\n",sep = "")
		#cat("\t","Incomplete Block Design","\n",sep = "") 
		cat("\t","Balanced Incomplete Block Design","\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Treatments = ", length(generate[[1]]), "\n",sep = "")
		cat("\t","Plots per Block (Block Size) = ", blkSize, "\n",sep = "")
		cat("\t","Number of Blocks = ", numBlks, "\n",sep = "")
		cat("\t","Number of Replicates = ", r, "\n",sep = "")
		cat("\t","Lambda = ", lambda, "\n\n",sep = "")
		#cat("Results of Randomization:\n")
		#printDataFrame(randomize)
		cat("\t","Number of Field Row = ", numFieldRow, "\n",sep = "")
		cat("\t","Number of Field Column = ", ncol(plan[[1]]), "\n",sep = "")
	}

	if (!is.null(file)) {
		tempFile <- strsplit(file, split = "\\.")[[1]]
		tempExt <- tolower(tempFile[length(tempFile)])
		if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
		newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
		newFile <- paste(newFile, tempExt, sep = ".")
		if (tempExt == "csv") { write.csv(randomize, file = newFile, row.names = FALSE)
		} else { save(randomize, file = newFile) }
	} else {
		cat("Results of Randomization:\n")
		printDataFrame(randomize)
	}
	return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum)))
}

