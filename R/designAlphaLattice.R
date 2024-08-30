# -----------------------------------------------------------------------------------------
# designAlphaLattice: Generate randomization for alpha lattice design using the 4 series
#             formulated by Patterson and Williams
# Created by: Alaine A. Gulles 07.16.2013 for International Rice Research Institute
# Modified by: Alaine A. Gulles 08.29.2014
# Note: Uses the DiGGer function from package DiGGer
# -----------------------------------------------------------------------------------------

designAlphaLattice <- function(generate, blksize, r = 2, trial = 1, rowPerBlk, rowPerRep, numFieldRow,
                               serpentine = FALSE, file = NULL) UseMethod("designAlphaLattice")

designAlphaLattice.default <- function(generate, blksize, r = 2, trial = 1, rowPerBlk, rowPerRep, numFieldRow, serpentine = FALSE, file = NULL) {
     
     # --- check entry --- #
     
     # check if r greater than 2    
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.")}
     if (length(generate[[1]]) == 1) { tempComb <- FactorList(generate) } else { tempComb <- generate }
     
     if (rowPerBlk == 1) serpentine <- FALSE
     
     # determine the total number of experimental units
     if ((r * length(tempComb[[1]])) > 1500) { stop("The maximum number of experimental units that can be generated is 1500.") }
     
     # check if the number of treatment is divisible by the block size
     if (length(tempComb[[1]])%%blksize != 0) { stop("The number of treatments should be divisible by the number of plots per block.") }
     numBlk <- length(tempComb[[1]])/blksize                # determine the number of block per replicate
     
     if (blksize%%rowPerBlk != 0) { stop("The number of plots per block should be divisible by the number of rows per block.") }
     colPerBlk <- blksize/rowPerBlk                         # determine the number of columns per block with a replicate
     
     # check if # of rows per replicate is divisible by the # of rows per block
     if (rowPerRep%%rowPerBlk != 0) { stop("The number of rows per replicate should be divisible by the number of rows per block.") }
     
     # check if the quotient of # of rows per replicate and # of rows per block is a factor to number of blocks per replicate
     if (!((rowPerRep/rowPerBlk) %in% allFactors(numBlk))) { stop("The quotient of the number of rows in each replicate and number of rows in each block should be a factor of the number of blocks per replicate.") }
     
     # check if the treatment is divisible by the number of rows per Rep
     if (length(tempComb[[1]])%%rowPerRep != 0) { stop("The number of treatment should be divisible by the number of rows in each replicate.") }
     colPerRep <- length(tempComb[[1]])/rowPerRep           # determine the number of columns per replicate
     
          
     if (numFieldRow%%rowPerRep != 0) { stop("The total number of plots should be divisible by the number of field rows.") }
     numRepRow <- numFieldRow/rowPerRep                     # determine the number of rep along the length of the field layout
     
     if (!(numRepRow %in% allFactors(r))) { stop("The quotient of the number of field rows and number of rows in each replicate should be a factor of the number of the replicates.") }
     
     if((length(tempComb[[1]])*r)%%numFieldRow != 0) { stop("The total number of plots should be divisible by the number of field rows.") } 
     numFieldCol <- (length(tempComb[[1]])*r)/numFieldRow    # determine the number of field column in the experiment
          
     numRepCol <- numFieldCol/colPerRep                     # determine the number of blocks along the length of each replicate
     numBlkRow <- rowPerRep/rowPerBlk                       # determine the number of blocks along the length of each replicate
     numBlkCol <- colPerRep/colPerBlk                       
     
     randomize <- NULL
     plan <- list()
     plotNum <- NULL
     blkNum <- NULL

     if (r * length(tempComb[[1]]) != numFieldRow * numFieldCol) { stop("Total of plots cannot be accomodated in the field experiment.") }
     
     for (i in (1:trial)) {
          
          result <- try(temp <- DiGGer(NumberOfTreatments = length(tempComb[[1]]), 
                                       RowsInDesign = numFieldRow, 
                                       ColumnsInDesign = numFieldCol, 
                                       RowsInReplicate = rowPerRep, 
                                       ColumnsInReplicate = colPerRep,
                                       BlockIn2D = c(rowPerBlk,colPerBlk), 
                                       #RowsInBlockSequence = c(rowPerRep, rowPerBlk),
                                       #ColumnsInBlockSequence = c(colPerRep, colPerBlk),
                                       TreatmentName = tempComb[[1]]), silent = TRUE)
          
          
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          setCorrelation(temp, phasenumber = c(1, rowPerRep, colPerRep)) ### --- added by AAGulles c/o VIBartolome 26May2014
          
          capture.output(run(temp))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          result_mat <- getDesign(temp)
          #des.plot(result_mat, rstr = "", cstr = "", bdef = cbind(rowPerBlk, colPerBlk), bcol = 4)
          #des.plot(result_mat, rstr = "", cstr = "", bdef = cbind(rowPerRep, colPerRep), new = FALSE)          
          plan[[i]] <- matrix(print(temp, option = "list")$ID, nrow(result_mat), ncol(result_mat))
          
          if (i == 1) {
               tempfbook <- print(temp, option = "list")
               blkNumTemp <- matrix(0, rowPerRep, colPerRep)
               blkPlotNumTemp <- matrix(0, rowPerRep, colPerRep)
               
               
               #tempPlotNum <- matrix(1:length(tempComb[[1]]), rowPerRep, colPerRep, byrow = TRUE)
               #if (serpentine) { for (k in seq(2, rowPerRep, by = 2)) { tempPlotNum[k,] <- rev(tempPlotNum[k,]) }}
               
               tempPlotNumBlk <- matrix(1:blksize, rowPerBlk, colPerBlk, byrow = TRUE)
               if (serpentine) { for (k in seq(2, rowPerBlk, by = 2)) { tempPlotNumBlk[k,] <- rev(tempPlotNumBlk[k,]) }}
               
               blkcode <- 1
               index <- tempPlotNumBlk
               for (j in 1:numBlkRow) {
                    for (k in 1:numBlkCol) {
                         rowIndexLL <- (j * rowPerBlk) - rowPerBlk + 1
                         rowIndexUL <- rowIndexLL + rowPerBlk - 1
                         colIndexLL <- (k * colPerBlk) - colPerBlk + 1
                         colIndexUL <- colIndexLL + colPerBlk - 1
                         blkNumTemp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <-  blkcode
                         blkPlotNumTemp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- index
                         blkcode <- blkcode + 1
                         index <- index + blksize
                    }
               }
               blkNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               plotNum <- matrix(as.numeric(paste(tempfbook$REP,paste(c(rep(0, max(nchar(1:length(tempComb[[1]]))))), collapse = ""), sep = "")), nrow(plan[[i]]), ncol(plan[[i]]))
                                        
               for (j in 1:numRepRow) {
                    for (k in 1:numRepCol) {
                         rowIndexLL <- (j * rowPerRep) - rowPerRep + 1
                         rowIndexUL <- rowIndexLL + rowPerRep - 1
                         colIndexLL <- (k * colPerRep) - colPerRep + 1
                         colIndexUL <- colIndexLL + colPerRep - 1
                         #plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL]+tempPlotNum 
                         plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL]+blkPlotNumTemp 
                         blkNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- blkNumTemp
                    }
               }
               
          }
          
          tempFieldOrder <- merge(as.data.frame.table(blkNum), as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("Block","PlotNum"))
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"]) 
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          colnames(tempFieldOrder) <- gsub("Freq", "", colnames(tempFieldOrder))
          randomize <- rbind(randomize, cbind(Trial = i, merge(print(temp, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
          
     }

     dimnames(plotNum) <- dimnames(plan[[1]])
     names(plan) <- paste("Trial", 1:trial, sep = "")
     randomize <- randomize[,c("Trial", "REP", "Block", "ID", "PlotNum", "ROW", "RANGE")]
     names(randomize) <- c("Trial", "Rep", "Block", names(tempComb)[1], "PlotNum", "FieldRow", "FieldCol")
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum),]
     randomize[,"Trial"] <- factor(randomize[,"Trial"])
	randomize[,"Rep"] <- factor(randomize[,"Rep"])
	randomize[,"Block"] <- factor(randomize[,"Block"])
	randomize[,"FieldRow"] <- as.numeric(randomize[,"FieldRow"])
	randomize[,"FieldCol"] <- as.numeric(randomize[,"FieldCol"])
	rownames(randomize) <- 1:nrow(randomize)

	
	cat(toupper("Design Properties:"),"\n",sep = "")
	cat("\t","Incomplete Block Design","\n",sep = "") 
	cat("\t","Alpha Lattice Design","\n\n",sep = "") 
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
	cat("\t","Number of Replicates = ", r, "\n",sep = "")
	cat("\t","Number of Plots per Block = ", blksize, "\n",sep = "")
	cat("\t","Number of Blocks per Replicate = ", length(tempComb[[1]])/blksize, "\n\n",sep = "")
     
     cat("\t","Number of Field Rows = ", numFieldRow, "\n",sep = "")
     cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")
     
	#cat("Results of Randomization:\n")
	#printDataFrame(randomize)
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
	return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum, blkNum = blkNum)))
}

