# -------------------------------------------------------------------------------------
# designSplit: Generate randomization for split plot family design.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.17.2013
# -------------------------------------------------------------------------------------

designSplit <- function(main, sub, ssub = NULL, sssub = NULL, r = NULL, trial = 1, design = c('crd', 'rcbd', 'lsd'), 
                        numFieldRow = 1, rowPerBlk = 1, rowPerMain = 1, rowPerSub = 1, rowPerSubSub = 1, serpentine = FALSE, display = TRUE, file = NULL) UseMethod("designSplit")

designSplit.default <- function(main, sub, ssub = NULL, sssub = NULL, r = NULL, trial = 1, design = c('crd', 'rcbd', 'lsd'), 
                                numFieldRow = 1, rowPerBlk = 1, rowPerMain = 1, rowPerSub = 1, rowPerSubSub = 1, serpentine = FALSE, display = TRUE, file = NULL) {

	if (missing(main)) { stop("The argument 'main' is missing.") }
	if (missing(sub)) { stop("The argument 'sub' is missing.") }

	design <- match.arg(design)
	factorList <- c(main, sub, ssub, sssub)
	if (length(main[[1]]) == 1) { trmtLevel <- main[[1]]; factorList[[1]] <- 1:main[[1]] } else { trmtLevel <- length(main[[1]]) }
     
     ## if (numFieldRow == 1 || rowPerBlk == 1) serpentine <- FALSE
	if (numFieldRow == 1 && design != "lsd") serpentine <- FALSE

	# ---------------------------------------------------------------------
	# randomizing and layouting the main plot factor
	# ---------------------------------------------------------------------
	     
     if (design == "crd")  { 
          #rowPerBlk <- 1
          tempNumRow <- numFieldRow/(rowPerMain*rowPerSub*rowPerSubSub)
          tempNumCol <- (trmtLevel[1]*r)/tempNumRow
          capture.output(result <- designCRD(main, r, trial, numFieldRow = tempNumRow, display = FALSE))
          book <- result$fieldbook
     }
     
	if (design == "rcbd") { 
          numBlkRow <- numFieldRow/(rowPerBlk*rowPerMain*rowPerSub*rowPerSubSub)
          numBlkCol <- r/numBlkRow
          capture.output(result <- designRCBD(main, r, trial, numFieldRow = numFieldRow/(rowPerMain*rowPerSub*rowPerSubSub), rowPerBlk, display = FALSE))
          book <- result$fieldbook
	}
	if (design == "lsd")  { 
	     #rowPerBlk <- 1
          capture.output(result <- designLSD(main, trial, display = FALSE))		 
          book <- result$fieldbook
	}

     # ---------------------------------------------------------------------
	# randomizing and layouting the subplot factor within each main plot
	# ---------------------------------------------------------------------
     
	if (length(sub[[1]]) == 1) { trmtLevel <- c(trmtLevel, sub[[1]]); factorList[[2]] <- 1:sub[[1]] } else { trmtLevel <- c(trmtLevel, length(sub[[1]])) }
	colPerMain <- trmtLevel[2]/rowPerMain
     
	# for rcbd
     if (design == "rcbd") {
          rowPerBlk <- rowPerBlk*rowPerMain
          colPerBlk <- prod(trmtLevel)/rowPerBlk
          tempNumRow <- rowPerBlk * numBlkRow
          tempNumCol <- colPerBlk*numBlkCol
          totalPlots <- prod(trmtLevel)
     }
     if (design == "crd"){
        tempNumRow <- rowPerMain * tempNumRow
        tempNumCol <- colPerMain * tempNumCol
        totalPlots <- tempNumRow * tempNumCol
     }
	if (design == "lsd") {
	     tempNumRow <- rowPerMain * trmtLevel[1]
	     tempNumCol <- colPerMain * trmtLevel[1]
	     totalPlots <- prod(trmtLevel)*trmtLevel[1]
	     book$FieldRow <- book$Row
	     book$FieldCol <- book$Column
	} 
	
     randomize <- NULL
     plan <- list()
     plotNum <- NULL
	for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(GenerateFactor(sub)[,1], length(GenerateFactor(sub)[,1])))) }
	colnames(randomize)[ncol(randomize)] <- names(sub)[1]
	
	randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol),]
	randomize <- cbind(randomize, withinFCol = 1:colPerMain)
	randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol, randomize$withinFCol),]
	randomize <- cbind(randomize, withinFRow = 1:rowPerMain)
	randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$withinFRow, randomize$FieldCol, randomize$withinFCol),]
	randomize <- cbind(randomize, newFieldRow = rep(1:tempNumRow, each = tempNumCol))
	randomize <- randomize[order(randomize$Trial, randomize$FieldCol, randomize$withinFCol, randomize$newFieldRow),]
	randomize <- cbind(randomize, newFieldCol = rep(1:tempNumCol, each = tempNumRow))
     if (design == "crd" || design == "lsd") { randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),] 
     #} else {  randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$newFieldRow, randomize$newFieldCol),] }
     } else {  randomize <- randomize[order(randomize$Trial, randomize$Rep, randomize$newFieldRow, randomize$newFieldCol),] }
     randomize$tPNum <- 1:totalPlots
     if (design == "rcbd") { randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),] }
	randomize[,"trmtLabel"] <- paste(randomize[,names(main)[1]], randomize[,names(sub)[1]])
     
	for (i in (1:trial)) {
	     plan[[i]] <- matrix(randomize[randomize[,"Trial"] == 1, "trmtLabel"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
	}
     
	plotNum <- matrix(randomize[randomize[,"Trial"] == 1, "tPNum"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
     
	if (serpentine) {
          if (design == "rcbd") {
               if (rowPerBlk > 1){
                    for (i in (1:numBlkRow)) {
                         if (i == 1) { start <- 1; last <- rowPerBlk } else { start <- start + rowPerBlk; last <- last + rowPerBlk }
                         for (j in seq((start + 1), last, by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }
                    }
               }
          } else { for (j in seq(2, nrow(plotNum), by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }}
     }
	
     if (design == "rcbd") {
          #blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Block"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
          blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Rep"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
          blockNum <- matrix(as.numeric(paste(blockNum, paste(c(rep(0,max(nchar(plotNum)))), collapse = ""), sep = "")), nrow = tempNumRow, ncol = tempNumCol, byrow = FALSE)
          plotNum <- blockNum + plotNum
     }
     
     dimnames(plotNum) <- list(1:nrow(plotNum), 1:ncol(plotNum))
	randomize <- merge(randomize, as.data.frame.table(plotNum), by.x = c("newFieldRow", "newFieldCol"), by.y = c("Var1", "Var2"))
	randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
	
	if (design == "crd") {  
          randomize <- randomize[order(randomize$trmtLabel, randomize$Trial),]
          randomize$Rep <- 1:r
          randomize <- randomize[order(randomize$Trial, randomize$tPNum),]
	     randomize <- randomize[,c("trmtLabel", "Trial", "Rep",names(main)[[1]], names(sub)[1], "Freq", "newFieldRow", "newFieldCol")]
	}
     
     if (design == "rcbd") { randomize <- randomize[,c("trmtLabel", "Trial", "Rep",names(main)[[1]], names(sub)[1], "Freq", "newFieldRow", "newFieldCol")] }
     if (design == "lsd") { randomize <- randomize[,c("trmtLabel", "Trial", "Row", "Column", names(main)[[1]], names(sub)[1], "Freq", "newFieldRow", "newFieldCol")] }
     
     names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")	
     dimnames(plotNum) <- dimnames(plan[[1]])
     factorList[1] <- FactorList(main)
	factorList[2] <- FactorList(sub)
     
	# -----------------------------------------------------------------------------
	# randomizing and layouting the sub-subplot factor within each sub plot
	# -----------------------------------------------------------------------------
	
	if (!is.null(ssub)) {

	     if (length(ssub[[1]]) == 1) { trmtLevel <- c(trmtLevel, ssub[[1]]); factorList[[3]] <- 1:ssub[[1]] } else { trmtLevel <- c(trmtLevel, length(ssub[[1]])) }
	     colPerSub <- trmtLevel[3]/rowPerSub
          
          if (design == "crd" || design == "lsd") {
               tempNumRow <- rowPerSub * tempNumRow
               tempNumCol <- colPerSub * tempNumCol
               totalPlots <- tempNumRow * tempNumCol
          } else {
	          rowPerBlk <- rowPerBlk*rowPerSub
	          colPerBlk <- prod(trmtLevel)/rowPerBlk
               tempNumRow <- rowPerSub * tempNumRow
               tempNumCol <- colPerSub * tempNumCol
               totalPlots <- prod(trmtLevel)
	     }
	               
          book <- randomize
		randomize <- NULL
		plan <- list()
		plotNum <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(GenerateFactor(ssub)[,1], length(GenerateFactor(ssub)[,1])))) }
		colnames(randomize)[ncol(randomize)] <- names(ssub)[1]
		
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol),]
		randomize <- cbind(randomize, withinFCol = 1:colPerSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, withinFRow = 1:rowPerSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$withinFRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, newFieldRow = rep(1:tempNumRow, each = tempNumCol))
          randomize <- randomize[order(randomize$Trial, randomize$FieldCol, randomize$withinFCol, randomize$newFieldRow),]
          randomize <- cbind(randomize, newFieldCol = rep(1:tempNumCol, each = tempNumRow))
          if (design == "crd" || design == "lsd") { randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]     
          } else { randomize <- randomize[order(randomize$Trial, randomize$Rep, randomize$newFieldRow, randomize$newFieldCol),] }
		randomize$tPNum <- 1:totalPlots
		if (design == "rcbd") randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
		randomize[,"trmtLabel"] <- paste(randomize[,"trmtLabel"], randomize[,names(ssub)[1]])
		
		for (i in (1:trial)) {
		     plan[[i]] <- matrix(randomize[randomize[,"Trial"] == 1, "trmtLabel"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		     dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
		}
          
		plotNum <- matrix(randomize[randomize[,"Trial"] == 1, "tPNum"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		
		if (serpentine) {
               if (design == "rcbd") {
                    if (rowPerBlk > 1) {
                         for (i in (1:numBlkRow)) {
                              if (i == 1) { start <- 1; last <- rowPerBlk } else { start <- start + rowPerBlk; last <- last + rowPerBlk }
                              for (j in seq((start + 1), last, by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }
                         }     
                    }
               } else { for (j in seq(2, nrow(plotNum), by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }}
		}
		
		if (design == "rcbd") {
		     blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Rep"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		     blockNum <- matrix(as.numeric(paste(blockNum, paste(c(rep(0,max(nchar(plotNum)))), collapse = ""), sep = "")), nrow = tempNumRow, ncol = tempNumCol, byrow = FALSE)
		     plotNum <- blockNum + plotNum
		}
		dimnames(plotNum) <- list(1:nrow(plotNum), 1:ncol(plotNum))
		randomize <- merge(randomize, as.data.frame.table(plotNum), by.x = c("newFieldRow", "newFieldCol"), by.y = c("Var1", "Var2"))
		randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
	     if (design == "crd") { 
               randomize <- randomize[order(randomize$trmtLabel, randomize$Trial),]
               randomize$Rep <- 1:r
               randomize <- randomize[order(randomize$Trial, randomize$tPNum),]
               randomize <- randomize[,c("trmtLabel", "Trial", "Rep", names(main)[[1]], names(sub)[1], names(ssub)[1], "Freq", "newFieldRow", "newFieldCol")] 
	     }
		if (design == "rcbd") { randomize <- randomize[,c("trmtLabel", "Trial", "Rep",names(main)[[1]], names(sub)[1], names(ssub)[1], "Freq", "newFieldRow", "newFieldCol")] }
	     if (design == "lsd") { randomize <- randomize[,c("trmtLabel", "Trial", "Row", "Column", names(main)[[1]], names(sub)[1], names(ssub)[1], "Freq", "newFieldRow", "newFieldCol")] }
		names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
		dimnames(plotNum) <- dimnames(plan[[1]])
		
	     factorList[3] <- FactorList(ssub)
	}
     
	# -----------------------------------------------------------------------------
	# randomizing and layouting the sub-sub-subplot factor within each sub-subplot
	# -----------------------------------------------------------------------------
	     
	if (!is.null(sssub)) {
		if (length(sssub[[1]]) == 1) { trmtLevel <- c(trmtLevel, sssub[[1]]); factorList[[4]] <- 1:sssub[[1]] } else { trmtLevel <- c(trmtLevel, length(sssub[[1]])) }
		colPerSubSub <- trmtLevel[4]/rowPerSubSub
          
		if (design == "crd" || design == "lsd") {
		     tempNumRow <- rowPerSubSub * tempNumRow
		     tempNumCol <- colPerSubSub * tempNumCol
		     totalPlots <- tempNumRow * tempNumCol
		} else {
		     rowPerBlk <- rowPerBlk*rowPerSubSub
		     colPerBlk <- prod(trmtLevel)/rowPerBlk
		     tempNumRow <- rowPerSubSub * tempNumRow
		     tempNumCol <- colPerSubSub * tempNumCol
		     totalPlots <- prod(trmtLevel)
		}
		
		book <- randomize
		randomize <- NULL
		plan <- list()
		plotNum <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(GenerateFactor(sssub)[,1], length(GenerateFactor(sssub)[,1])))) }
		colnames(randomize)[ncol(randomize)] <- names(sssub)[1]
		
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol),]
		randomize <- cbind(randomize, withinFCol = 1:colPerSubSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, withinFRow = 1:rowPerSubSub)
		randomize <- randomize[order(randomize$Trial, randomize$FieldRow, randomize$withinFRow, randomize$FieldCol, randomize$withinFCol),]
		randomize <- cbind(randomize, newFieldRow = rep(1:tempNumRow, each = tempNumCol))
		randomize <- randomize[order(randomize$Trial, randomize$FieldCol, randomize$withinFCol, randomize$newFieldRow),]
		randomize <- cbind(randomize, newFieldCol = rep(1:tempNumCol, each = tempNumRow))
          if (design == "crd" || design == "lsd") {  randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
          } else { randomize <- randomize[order(randomize$Trial, randomize$Rep, randomize$newFieldRow, randomize$newFieldCol),] }
		randomize$tPNum <- 1:totalPlots
          if (design == "rcbd") { randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),] }
		randomize[,"trmtLabel"] <- paste(randomize[,"trmtLabel"], randomize[,names(sssub)[1]])
		
		for (i in (1:trial)) {
		     plan[[i]] <- matrix(randomize[randomize[,"Trial"] == 1, "trmtLabel"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		     dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
		}
		
		plotNum <- matrix(randomize[randomize[,"Trial"] == 1, "tPNum"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		
		if (serpentine) {
               if (design == "rcbd") {
                    if (rowPerBlk > 1) {
                         for (i in (1:numBlkRow)) {
                              if (i == 1) { start <- 1; last <- rowPerBlk } else { start <- start + rowPerBlk; last <- last + rowPerBlk }
                              for (j in seq((start + 1), last, by = 2)) { plotNum[j, ] <- rev(plotNum[j, ]) }
                         }
                    }
               } else { for (j in seq(2, nrow(plotNum), by = 2)) { plotNum[j,] <- rev(plotNum[j,]) }}
		}
		
		if (design == "rcbd") {
		     blockNum <- matrix(randomize[randomize[,"Trial"] == 1, "Rep"], nrow = tempNumRow, ncol = tempNumCol, byrow = TRUE)
		     blockNum <- matrix(as.numeric(paste(blockNum, paste(c(rep(0,max(nchar(plotNum)))), collapse = ""), sep = "")), nrow = tempNumRow, ncol = tempNumCol, byrow = FALSE)
		     plotNum <- blockNum + plotNum
		}
		dimnames(plotNum) <- list(1:nrow(plotNum), 1:ncol(plotNum))
		randomize <- merge(randomize, as.data.frame.table(plotNum), by.x = c("newFieldRow", "newFieldCol"), by.y = c("Var1", "Var2"))
		randomize <- randomize[order(randomize$Trial, randomize$newFieldRow, randomize$newFieldCol),]
          
          if (design == "crd") {
               randomize <- randomize[order(randomize$trmtLabel, randomize$Trial),]
               randomize$Rep <- 1:r
               randomize <- randomize[order(randomize$Trial, randomize$tPNum),]
               randomize <- randomize[,c("trmtLabel", "Trial", "Rep", names(main)[[1]], names(sub)[1], names(ssub)[1], names(sssub)[1], "Freq", "newFieldRow", "newFieldCol")] 
          }
		if (design == "rcbd") { randomize <- randomize[,c("trmtLabel", "Trial", "Rep",names(main)[[1]], names(sub)[1], names(ssub)[1], names(sssub)[1], "tPNum", "newFieldRow", "newFieldCol")] }
		if (design == "lsd") { randomize <- randomize[,c("trmtLabel", "Trial", "Row", "Column", names(main)[[1]], names(sub)[1], names(ssub)[1], names(sssub)[1], "tPNum", "newFieldRow", "newFieldCol")] }
     	names(randomize)[(ncol(randomize) - 2):ncol(randomize)] <- c("PlotNum", "FieldRow", "FieldCol")
		dimnames(plotNum) <- dimnames(plan[[1]])
		factorList[4] <- FactorList(sssub)
	}
     
	numFactor <- length(trmtLevel)
	randomize <- randomize[,-I(c(match("trmtLabel", names(randomize))))]
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum),]
	rownames(randomize) <- 1:nrow(randomize)
     
     # to include mainplot, subplot, ssubplot and sssubplot column
     Mainplot <- randomize[,names(factorList)[1]]
	Subplot <- randomize[,names(factorList)[2]]
     if (length(factorList) > 2) {  
          SubSubplot <- randomize[,names(factorList)[3]]
          if (length(factorList) > 3) { SubSubSubplot <- randomize[,names(factorList)[4]]  } else { SubSubSubplot <- NULL }  
     } else { 
          SubSubplot <- NULL
          SubSubSubplot <- NULL 
     }
     
     if (is.null(SubSubplot) && is.null(SubSubSubplot)) {
          randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
                             Mainplot = Mainplot, Subplot = Subplot, randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
     }          

	if (!is.null(SubSubplot) && is.null(SubSubSubplot)) {
	     randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
	                        Mainplot = Mainplot, Subplot = Subplot, SubSubplot = SubSubplot, randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
	}          
	
	if (!is.null(SubSubplot) && !is.null(SubSubSubplot)) {
	     randomize <- cbind(randomize[c(1:(match(names(factorList)[1], names(randomize)))-1)],
	                        Mainplot = Mainplot, Subplot = Subplot, SubSubplot = SubSubplot, SubSubSubplot = SubSubSubplot, 
                             randomize[match(names(factorList)[1], names(randomize)):ncol(randomize)])
	}          
	
     
	if (display) {
		splitLabel <- c("Split Plot Design", "Split-Split Plot Design", "Split-Split-Split Plot Design")
		designLabel <- c("Completely Randomized Design", "Randomized Complete Block Design", "Latin Square Design")
		factorLabel <- c("Mainplot", "Subplot", "Sub-subplot", "Sub-sub-subplot")

		cat(toupper("Design Properties:"),"\n",sep = "")
		cat("\t",splitLabel[numFactor-1],"\n",sep = "") 
		cat("\t",designLabel[match(design, c("crd", "rcbd", "lsd"))],"\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		if (design == "crd") cat("\t","Number of Replicates = ", r, "\n",sep = "")
		if (design == "rcbd") cat("\t","Number of Replicates = ", r, "\n",sep = "")
		for (i in (1:length(factorList))) {
			cat("\t",factorLabel[i]," Factor = ", names(factorList)[i], "\n",sep = "")
			cat("\t","Levels = ", sep = "")
			if (trmtLevel[i] <= 5) { cat(paste(factorList[[i]], collapse = ", ", sep = ""), sep = "")
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