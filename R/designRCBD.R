# -------------------------------------------------------------------------------------
# designRCBD: Generate randomization for randomized complete block design (RCBD)
#             for single factor or factorial experiments.
# Created by: Alaine A. Gulles 04.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.11.2012
# -------------------------------------------------------------------------------------

designRCBD <- function(generate, r = 2, trial = 1, numFieldRow = 1, rowPerBlk = 1, serpentine = FALSE, display = TRUE, file = NULL) UseMethod("designRCBD")

designRCBD.default <- function(generate, r = 2, trial = 1, numFieldRow = 1, rowPerBlk = 1, serpentine = FALSE, display = TRUE, file = NULL) {

	if (is.null(trial) || trial < 1 || is.character(trial) || length(trial) > 1) { stop("The argument 'trial' should be a single value greater than or equal to 1.") }
	if (is.null(r) || r < 2 || is.character(r) || length(r) > 1) { stop("The argument 'r' should be a single value greater than or equal to 2.") }
	if (missing(generate)) { stop("The argument 'generate' is missing.") }
	if (!is.list(generate)) { stop("The argument 'generate' must be a list.") }
     if (rowPerBlk > numFieldRow) { stop("Number of field row should be equal to greater than the number or row per rep.") }
	
	tempComb <- GenerateFactor(generate, times = 1)
	randomize <- NULL
     #plotNum <- list()
     plan <- list()
     
     if (rowPerBlk == 1) serpentine <- FALSE
     
     numBlkRow <- numFieldRow/rowPerBlk
	numBlkCol <- r/numBlkRow
	colPerBlk <- (nrow(tempComb)/rowPerBlk)
     
     if (nrow(tempComb)%%rowPerBlk != 0) { stop("Total number of plots per replicate should be divisible by the number of rows within replicate.") }
	if ((nrow(tempComb)*r)%%numFieldRow != 0) { stop("Total number of plots should be divisible by the number of field rows.") }

	for (i in (1:trial)) {
          plan[[i]] <- matrix(0, nrow = numFieldRow, ncol = (nrow(tempComb)*r)/numFieldRow)
          if (i == 1) { plotNum <- plan[[i]] }
		for (j in (1:r)) {
		     tempPlan <- NULL
		     tempPlotNum <- NULL
			temp <- data.frame(Trial = as.character(i), Rep = as.character(j), tempComb, tempPlotNum = sample(nrow(tempComb), nrow(tempComb), replace = FALSE))
			temp <- temp[order(temp[,"tempPlotNum"]),]
			plotLabel <- as.numeric(paste(j, paste(c(rep(0, max(nchar(1:nrow(tempComb))))), collapse = ""), sep = ""))+1:nrow(tempComb)
			if (ncol(tempComb) > 1) { trmtLabel <- eval(parse(text = paste("paste(temp[,'", paste(names(temp)[3:(ncol(temp)-1)], collapse = "'],' ',temp[,'", sep = ""),"'], sep = '')", sep = "")))
			} else { trmtLabel <- temp[,3] }
			tempPlan <- matrix(trmtLabel, nrow = rowPerBlk, ncol = nrow(tempComb)/rowPerBlk, byrow = TRUE)
		     tempPlotNum <- matrix(plotLabel, nrow = rowPerBlk, ncol = nrow(tempComb)/rowPerBlk, byrow = TRUE)
			if (serpentine) { for (k in seq(2, rowPerBlk, by = 2)) { tempPlotNum[k, ] <- rev(tempPlotNum[k, ]) }}
     	     if (j%%numBlkCol != 0) { colIndex <- j%%numBlkCol } else { colIndex <- numBlkCol }
		     if (colIndex == 1) { colIndexLL <- colIndex  } else { colIndexLL <- colIndexLL + colPerBlk }
		     colIndexUL <- colIndexLL + colPerBlk - 1
		     rowIndex <- ceiling(j/numBlkCol)
		     rowIndexLL <- (rowIndex * rowPerBlk) - rowPerBlk + 1
		     rowIndexUL <- rowIndexLL + rowPerBlk - 1 
               plan[[i]][rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlan
               if (i == 1) plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlotNum
		     tempFieldOrder <- merge(as.data.frame.table(matrix(rowIndexLL:rowIndexUL,nrow = rowPerBlk, ncol = nrow(tempComb)/rowPerBlk, byrow = FALSE)), 
                                       as.data.frame.table(matrix(colIndexLL:colIndexUL,nrow = rowPerBlk, ncol = nrow(tempComb)/rowPerBlk, byrow = TRUE)), 
                                       by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
		     tempOrder <- merge(as.data.frame.table(tempPlan), as.data.frame.table(tempPlotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
		     randomize <- rbind(randomize, merge(cbind(temp, trmtLabel), merge(tempOrder, tempFieldOrder, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))[,3:6], by.x = "trmtLabel", by.y = "Freq.x.x"))
		}
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
          if (i == 1) dimnames(plotNum) <- dimnames(plan[[i]])
          
	}
     names(plan) <- paste("Trial", 1:trial, sep = "")
     randomize <- randomize[,-I(c(match(c("trmtLabel", "tempPlotNum"), names(randomize))))]
     names(randomize)[(ncol(randomize)-2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
	#randomize <- randomize[order(randomize$Trial, randomize$Rep, randomize$FieldRow, randomize$FieldCol),]
	randomize <- randomize[order(randomize$Trial, randomize$Rep, randomize$PlotNum),]
	rownames(randomize) <- 1:nrow(randomize)

	if (display) {
		cat(toupper("Design Properties:"),"\n",sep = "")
		if (ncol(tempComb) == 1) { cat("\t","Single Factor","\n",sep = "") } else { cat("\t","Factorial Design","\n",sep = "") }
		cat("\t","Randomized Complete Block Design","\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Replicates = ", r, "\n",sep = "")
		if (ncol(tempComb) == 1) {
			cat("\t","Treatment Name = ", names(tempComb)[1], "\n",sep = "")
			cat("\t","Treatment Levels = ", sep = "")
			if (nlevels(tempComb[,1]) <= 5) { cat(paste(levels(tempComb[,1]), collapse = ", ", sep = ""), sep = "")
			} else {
				cat(paste(levels(tempComb[,1])[1:3], collapse = ", ", sep = ""), sep = "")
				cat(paste(", ...,", levels(tempComb[,1])[nlevels(tempComb[,1])]), sep = "")
			}
			cat("\n\n")
		} else {
			for (i in (1:ncol(tempComb))) {
				cat("\t","Factor ",i," = ", names(tempComb)[i], "\n",sep = "")
				cat("\t","Levels = ", sep = "")
				if (nlevels(tempComb[,i]) <= 5) { cat(paste(levels(tempComb[,i]), collapse = ", ", sep = ""), sep = "")
				} else {
					cat(paste(levels(tempComb[,i])[1:3], collapse = ", ", sep = ""), sep = "")
					cat(paste(", ...,", levels(tempComb[,i])[nlevels(tempComb[,i])]), sep = "")
				}
				cat("\n")
			}
			cat("\n")
		}
		cat("\t","Number of Field Row = ", numFieldRow, "\n",sep = "")
		cat("\t","Number of Field Column = ", ncol(plan[[1]]), "\n",sep = "")
		#cat("Results of Randomization:\n")
		#printDataFrame(randomize)
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