# -----------------------------------------------------------------------------------------
# designRowColumn: Generate randomization for Row-Column Design
# Created by: Alaine A. Gulles 10.08.2013 for International Rice Research Institute 10.08.2013
# Modified by: Alaine A. Gulles 05.29.2013
# Note: Uses the DiGGer function from package DiGGer
# -----------------------------------------------------------------------------------------

designRowColumn <- function(generate, r = 2, trial = 1, rowPerRep, numFieldRow, serpentine = FALSE, file = NULL) UseMethod("designRowColumn")

designRowColumn.default <- function(generate, r = 2, trial = 1, rowPerRep, numFieldRow, serpentine = FALSE, file = NULL) {

     
     # --- CHECK INPUT
          
     # check if r greater than 2    
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.")}
     
     # create the levels of the treatment
     if (length(generate[[1]]) == 1) { tempComb <- FactorList(generate) } else { tempComb <- generate }
     
     # check if rowPerRep is greater than 1 or
     if (rowPerRep <= 1 || rowPerRep >= length(tempComb[[1]])) { stop("Number of row per replicate should not be equal to 1 or the number of treatments.") }
     
     # determine the total number of experimental units limitation of DiGGer package
     if ((r * length(tempComb[[1]])) > 1500) { stop("The maximum number of experimental units that can be generated is 1500.") }
     
     if ((length(tempComb[[1]])*r)%%numFieldRow != 0) { stop("The total number of plots should be divisible by the number of field rows.") }
     numFieldCol <- (length(tempComb[[1]])*r)/numFieldRow    # determine the number of field column in the experiment
     
     if (numFieldRow%%rowPerRep != 0) { stop("The total number of plots should be divisible by the number of field rows.") }
     numRepRow <- numFieldRow/rowPerRep                     # determine the number of rep along the length of the field layout
     
     if (!(numRepRow %in% allFactors(r))) { stop("The quotient of the number of field rows and number of rows in each replicate should be a factor of the number of the replicates.") }

     
     # determine the number of columns in each replicate
     if ((length(tempComb[[1]]) %% rowPerRep) != 0) { stop("The total number of treatment levels should be divisible by the number of replicate.") }
     colPerRep <- length(tempComb[[1]])/rowPerRep   # determine the number of columns per replicate
     
     if (numFieldCol%%colPerRep != 0) { stop("The total number of plots should be divisible by the number of field rows.") }
     numRepCol <- numFieldCol/colPerRep

     randomize <- NULL
     plan <- list()
     plotNum <- NULL
    

     if (r * length(tempComb[[1]]) != numFieldRow * numFieldCol) { stop("Total of plots cannot be accomodated in the field experiment.") }
     
     
     for (i in (1:trial)) {
          
          result <- try(temp <- DiGGer(NumberOfTreatments = length(tempComb[[1]]), 
                                        RowsInDesign = numFieldRow, 
                                        ColumnsInDesign = numFieldCol, 
                                        RowsInReplicate = rowPerRep, 
                                        ColumnsInReplicate = colPerRep,
                                        TreatmentName = tempComb[[1]]), silent = TRUE)
         
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               cat("Error in DiGGer:", msg, "\n")
          }
          
          setCorrelation(temp, phasenumber = c(1, rowPerRep, colPerRep)) ### --- added by AAGulles c/o VIBartolome 27May2014
                    
          capture.output(run(temp))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               cat("Error in DiGGer:", msg, "\n")
          }
          result_mat <- getDesign(temp)
          #des.plot(result_mat, rstr = "", cstr = "", bdef = cbind(rowPerRep, colPerRep))          
          plan[[i]] <- matrix(print(temp, option = "list")$ID, nrow(result_mat), ncol(result_mat))
          
          if (i == 1) {
               tempfbook <- print(temp, option = "list")
               plotNum <- matrix(as.numeric(paste(tempfbook$REP,paste(c(rep(0, max(nchar(1:length(tempComb[[1]]))))), collapse = ""), sep = "")), nrow(plan[[i]]), ncol(plan[[i]]))
               tempPlotNum <- matrix(1:length(tempComb[[1]]), rowPerRep, colPerRep, byrow = TRUE)
               if (serpentine) { for (k in seq(2, rowPerRep, by = 2)) { tempPlotNum[k,] <- rev(tempPlotNum[k,]) }}
               for (j in 1:numRepRow) {
                    for (k in 1:numRepCol) {
                         rowIndexLL <- (j * rowPerRep) - rowPerRep + 1
                         rowIndexUL <- rowIndexLL + rowPerRep - 1
                         colIndexLL <- (k * colPerRep) - colPerRep + 1
                         colIndexUL <- colIndexLL + colPerRep - 1
                         plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL]+tempPlotNum 
                    }
               }
          } ## end stmt -- if (i == 1)
          
          tempFieldOrder <- as.data.frame.table(plotNum)
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"]) 
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          colnames(tempFieldOrder)[3] <- "PlotNum"
          randomize <- rbind(randomize, cbind(Trial = i, merge(print(temp, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
     } ## end stmt -- for (i in (1:trial))

     dimnames(plotNum) <- dimnames(plan[[1]])
     names(plan) <- paste("Trial", 1:trial, sep = "")
     randomize$RowBlk <- randomize$ROW%%rowPerRep
     randomize[randomize[,"RowBlk"] == 0, "RowBlk"] <- rowPerRep
     randomize$ColBlk <- randomize$RANGE
     randomize <- randomize[,c("Trial", "REP", "RowBlk", "ColBlk", "ID", "PlotNum", "ROW", "RANGE")]
     names(randomize) <- c("Trial", "Rep", "RowBlk", "ColBlk", names(tempComb)[1], "PlotNum", "FieldRow", "FieldCol")
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum),]
     randomize[,"Trial"] <- factor(randomize[,"Trial"])
	randomize[,"Rep"] <- factor(randomize[,"Rep"])
	randomize[,"FieldRow"] <- as.numeric(randomize[,"FieldRow"])
	randomize[,"FieldCol"] <- as.numeric(randomize[,"FieldCol"])
	rownames(randomize) <- 1:nrow(randomize)

	
	cat(toupper("Design Properties:"),"\n",sep = "")
	cat("\t","Incomplete Block Design","\n",sep = "") 
	cat("\t","Row-Column Design","\n\n",sep = "") 
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
	cat("\t","Number of Replicates = ", r, "\n",sep = "")
     cat("\t","Number of Rows per Replicate = ", rowPerRep, "\n",sep = "")
     cat("\t","Number of Columns per Replicate = ", colPerRep, "\n",sep = "")
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
	return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum)))
}

