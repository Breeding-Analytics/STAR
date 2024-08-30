# -------------------------------------------------------------------------------------
# designStrip: Generate randomization for strip plot family design.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.13.2012
# -------------------------------------------------------------------------------------

designStrip <- function(vertical, horizontal, sub = NULL, ssub = NULL, r = 2, trial = 1, numFieldRow, 
                        rowPerMain = 1, rowPerSub = 1, serpentine = FALSE, display = TRUE, file = NULL) UseMethod("designStrip")

designStrip.default <- function(vertical, horizontal, sub = NULL, ssub = NULL, r = 2, trial = 1, numFieldRow, 
                                rowPerMain = 1, rowPerSub = 1, serpentine = FALSE, display = TRUE, file = NULL) {

     if (missing(vertical)) { stop("The argument 'vertical' is missing.") }
     if (missing(horizontal)) { stop("The argument 'horizontal' is missing.") }
     
     factorList <- c(FactorList(vertical),FactorList(horizontal))
     numItems <- prod(length(factorList[[1]]), length(factorList[[2]]))
     rowPerBlk <- length(factorList[[2]])
     colPerBlk <- length(factorList[[1]])
     if (numFieldRow%%(rowPerBlk*rowPerMain*rowPerSub) != 0) { stop("Total number of plots is not divisible by the number of field rows.") }
     numBlkRow <- numFieldRow/(rowPerBlk*rowPerMain*rowPerSub)
     numBlkCol <- r/numBlkRow
     
     randomize <- NULL
     plan <- list()
     plotNum <- list()
     
     # check for strip plot only
     for (i in (1:trial)) {
          plan[[i]] <- matrix(0, nrow = (numFieldRow/(rowPerMain*rowPerSub)), ncol = colPerBlk * numBlkCol)
          if (i == 1) { plotNum <- plan[[i]] }
          for (j in (1:r)) {
               vertical <- sample(factorList[[1]], length(factorList[[1]]))
               horizontal <- sample(factorList[[2]], length(factorList[[2]]))
               tempFieldBook <- data.frame(row.names = NULL, Trial = i, Block = j, rep(vertical, each = length(factorList[[2]])), rep(horizontal, length(factorList[[1]])))
               names(tempFieldBook)[3:4] <- c("Vertical", "Horizontal") 
               plotLabel <- as.numeric(paste(j, paste(c(rep(0, max(nchar(1:nrow(tempFieldBook))))), collapse = ""), sep = ""))+1:nrow(tempFieldBook)
               trmtLabel <- eval(parse(text = paste("paste(tempFieldBook[,'", paste(names(tempFieldBook)[3:4], collapse = "'],' ',tempFieldBook[,'", sep = ""),"'], sep = '')", sep = "")))
               tempPlan <- matrix(trmtLabel, nrow = rowPerBlk, ncol = colPerBlk, byrow = FALSE)
               tempPlotNum <- matrix(plotLabel, nrow = rowPerBlk, ncol = colPerBlk, byrow = TRUE)
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
               randomize <- rbind(randomize, merge(cbind(tempFieldBook, trmtLabel), merge(tempOrder, tempFieldOrder, by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))[,3:6], by.x = "trmtLabel", by.y = "Freq.x.x"))
          }
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
          if (i == 1) dimnames(plotNum) <- dimnames(plan[[i]])
     }
     names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")

	if (!is.null(sub)) {
		factorList <- c(factorList, FactorList(sub))
		numItems <- prod(numItems, length(factorList[[3]]))
          rowPerBlk <- rowPerBlk*rowPerMain
          colPerBlk <- numItems/rowPerBlk
          colPerMain <- length(factorList[[3]])/rowPerMain
		book <- randomize
          tempPlan <- plan
          tempPlotNum <- plotNum
		randomize <- NULL
          plan <- list()
          plotNum <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(factorList[[3]], length(factorList[[3]])))) }
          names(randomize)[ncol(randomize)] <- "sub"
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol),]
		randomize <- cbind(randomize, withinFCol = 1:colPerMain)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, withinFRow = 1:rowPerMain)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$withinFRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, newFieldRow = rep(1:(rowPerBlk*numBlkRow), each = colPerBlk*numBlkCol))
		randomize <- randomize[order(randomize$Trial, randomize$FieldCol, randomize$withinFCol, randomize$newFieldRow),]
		randomize <- cbind(randomize, newFieldCol = rep(1:(colPerBlk*numBlkCol), each = rowPerBlk*numBlkRow))
		randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$newFieldRow, randomize$newFieldCol),]
		randomize$tPNum <- 1:numItems
		randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
		randomize[,"trmtLabel"] <- paste(randomize[,"trmtLabel"], randomize[,"sub"])
          for (i in (1:trial)) {
               plan[[i]] <- matrix(randomize[randomize[,"Trial"] == 1, "trmtLabel"],
                                   nrow = rowPerBlk*numBlkRow,
                                   ncol = colPerBlk*numBlkCol, byrow = TRUE)
               dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
          }
		plotNum <- matrix(randomize[randomize[,"Trial"] == 1, "tPNum"],
		                    nrow = rowPerBlk*numBlkRow,
		                    ncol = colPerBlk*numBlkCol, byrow = TRUE)
          if (serpentine) {
               for (i in (1:numBlkRow)) {
                    if (i == 1) { start <- 1; last <- rowPerBlk } else { start <- start + rowPerBlk; last <- last + rowPerBlk }
                    for (j in seq((start + 1), last, by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }
               }
          }
		
		blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Block"],
		                   nrow = rowPerBlk*numBlkRow,
		                   ncol = colPerBlk*numBlkCol, byrow = TRUE)
          
		blockNum <- matrix(as.numeric(paste(blockNum, paste(c(rep(0,max(nchar(plotNum)))), collapse = ""), sep = "")),
		                   nrow = rowPerBlk*numBlkRow,
		                   ncol = colPerBlk*numBlkCol, byrow = FALSE)
          plotNum <- blockNum + plotNum
          dimnames(plotNum) <- list(1:nrow(plotNum), 1:ncol(plotNum))
		randomize <- merge(randomize, as.data.frame.table(plotNum), by.x = c("newFieldRow", "newFieldCol"), by.y = c("Var1", "Var2"))
          randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
		randomize <- randomize[,c("trmtLabel", "Trial", "Block","Vertical", "Horizontal", "sub", "tPNum", "newFieldRow", "newFieldCol")]
		names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
		dimnames(plotNum) <- dimnames(plan[[1]])
	}
     
     
	if (!is.null(ssub)) {
		factorList <- c(factorList, FactorList(ssub))
		numItems <- prod(numItems, length(factorList[[4]]))
		rowPerBlk <- rowPerBlk*rowPerSub
		colPerBlk <- numItems/rowPerBlk
		colPerSub <- length(factorList[[4]])/rowPerSub
		book <- randomize
		randomize <- NULL
		plan <- list()
		plotNum <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(factorList[[4]], length(factorList[[4]])))) }
		names(randomize)[ncol(randomize)] <- "ssub"
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol),]
          
		randomize <- cbind(randomize, withinFCol = 1:colPerSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, withinFRow = 1:rowPerSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$withinFRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, newFieldRow = rep(1:(rowPerBlk*numBlkRow), each = colPerBlk*numBlkCol))
		randomize <- randomize[order(randomize$Trial, randomize$FieldCol, randomize$withinFCol, randomize$newFieldRow),]
		randomize <- cbind(randomize, newFieldCol = rep(1:(colPerBlk*numBlkCol), each = rowPerBlk*numBlkRow))
		randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$newFieldRow, randomize$newFieldCol),]
		randomize$tPNum <- 1:numItems
		randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
		randomize[,"trmtLabel"] <- paste(randomize[,"trmtLabel"], randomize[,"ssub"])
		for (i in (1:trial)) {
		     plan[[i]] <- matrix(randomize[randomize[,"Trial"] == 1, "trmtLabel"],
		                         nrow = rowPerBlk*numBlkRow,
		                         ncol = colPerBlk*numBlkCol, byrow = TRUE)
		     dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
		}
		plotNum <- matrix(randomize[randomize[,"Trial"] == 1, "tPNum"],
		                  nrow = rowPerBlk*numBlkRow,
		                  ncol = colPerBlk*numBlkCol, byrow = TRUE)
		
          if (serpentine) {
               for (i in (1:numBlkRow)) {
                    if (i == 1) { start <- 1; last <- rowPerBlk } else { start <- start + rowPerBlk; last <- last + rowPerBlk }
                    for (j in seq((start + 1), last, by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }
               }
          }

          blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Block"],
		                   nrow = rowPerBlk*numBlkRow,
		                   ncol = colPerBlk*numBlkCol, byrow = TRUE)
		
		blockNum <- matrix(as.numeric(paste(blockNum, paste(c(rep(0,max(nchar(plotNum)))), collapse = ""), sep = "")),
		                   nrow = rowPerBlk*numBlkRow,
		                   ncol = colPerBlk*numBlkCol, byrow = FALSE)

          plotNum <- blockNum + plotNum
		dimnames(plotNum) <- list(1:nrow(plotNum), 1:ncol(plotNum))
		randomize <- merge(randomize, as.data.frame.table(plotNum), by.x = c("newFieldRow", "newFieldCol"), by.y = c("Var1", "Var2"))
		randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
		randomize <- randomize[,c("trmtLabel", "Trial", "Block","Vertical", "Horizontal", "sub", "ssub", "tPNum", "newFieldRow", "newFieldCol")]
		names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
		dimnames(plotNum) <- dimnames(plan[[1]])
	}
     randomize <- randomize[,-I(c(match("trmtLabel", names(randomize))))]
	colnames(randomize) <- c("Trial", "Block", names(factorList), "PlotNum", "FieldRow", "FieldCol")
     randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$FieldRow, randomize$FieldCol),]
     rownames(randomize) <- 1:nrow(randomize)
     
     # include vertical, horizontal, subplot and sub-subplot column in the fieldbook
     Vertical <- randomize[,names(factorList)[1]]
     Horizontal <- randomize[,names(factorList)[2]]
     if (length(factorList) > 2) {  
          Subplot <- randomize[,names(factorList)[3]]
          if (length(factorList) > 3) { SubSubplot <- randomize[,names(factorList)[4]]  } else { SubSubplot <- NULL }  
     } else { 
          Subplot <- NULL
          SubSubplot <- NULL 
     }
     
     if (is.null(Subplot) && is.null(Subplot)) {
          randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
                             Vertical = Vertical, Horizontal = Horizontal, randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
     }          
     
     if (!is.null(Subplot) && is.null(SubSubplot)) {
          randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
                             Vertical = Vertical, Horizontal = Horizontal, Subplot = Subplot, randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
     }          
     
     if (!is.null(Subplot) && !is.null(SubSubplot)) {
          randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
                             Vertical = Vertical, Horizontal = Horizontal, Subplot = Subplot, SubSubplot = SubSubplot, 
                             randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
     }          
     
     

 	if (display) {
		stripLabel <- c("Strip Plot Design", "Strip-Split Plot Design", "Strip-Split-Split Plot Design")
		factorLabel <- c("Vertical", "Horizontal", "Subplot", "Sub-subplot")
		cat(toupper("Design Properties:"),"\n",sep = "")
		cat("\t",stripLabel[length(factorList)-1],"\n\n",sep = "") 
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Blocks = ", r, "\n",sep = "")
		for (i in (1:length(factorList))) {
			cat("\t",factorLabel[i]," Factor = ", names(factorList)[i], "\n",sep = "")
			cat("\t","Levels = ", sep = "")
			if (length(factorList[[i]]) <= 5) { cat(paste(factorList[[i]], collapse = ", ", sep = ""), sep = "")
			} else {
				cat(paste(factorList[[i]][1:3], collapse = ", ", sep = ""), sep = "")
				cat(paste(", ...,", factorList[[i]][length(factorList[[i]])]), sep = "")
			}
			cat("\n")
		}
		cat("\n")
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