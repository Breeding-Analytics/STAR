# -------------------------------------------------------------------------
# designLatinizedAlpha: Generate randomization for latinized alpha design
# Created by: Alaine A. Gulles for International Rice Research Institute
# Note: Uses the DiGGer function from package DiGGer
# -------------------------------------------------------------------------

designLatinizedAlpha <- function(generate, blksize, r = 2, trial = 1, numFieldRow, serpentine = FALSE, file = NULL) UseMethod("designLatinizedAlpha")

designLatinizedAlpha.default <- function(generate, blksize, r = 2, trial = 1, numFieldRow, serpentine = FALSE, file = NULL) {
     
     # --- check entry --- #
     # check if number of replicate is greater than or equal to 2
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.")}
     
     # check the entry for the parameter generate
     if (length(generate[[1]]) == 1) { tempComb <- FactorList(generate) } else { tempComb <- generate }

     # determine the total number of experimental units
     if ((r * length(tempComb[[1]])) > 1500) { stop("The maximum number of experimental units that can be generated is 1500.") }
     
     # check if the number of treatment is divisible by the block size
     if (length(tempComb[[1]])%%blksize != 0) { stop("The number of treatments should be divisible by the number of plots per block.") }
     numBlk <- length(tempComb[[1]])/blksize                # determine the number of block per replicate
     
     # compute the number of plots per replicate
     numPlotsPerRep <- blksize * numBlk
     
     # determine if the number of plots per replicate is equal to the number of treatments 
     if (numPlotsPerRep != length(tempComb[[1]])) stop("ERROR: The number of plots per replicate should be equal to the number of treatment.")     
     
     # determine the case
     if (numFieldRow == numBlk) { case <- 2
     } else {
          if (numFieldRow == length(tempComb[[1]])) { case <- 1
          } else {
               if (numFieldRow > length(tempComb[[1]])) stop("ERROR: Number of field row should not be greater than the number of treatments.")
               if (numFieldRow%%blksize == 0) { case <- 1
               } else {
                    stop("ERROR: Number of field row should be divisible by the number of plots within a block.")     
               }
          }
     }
     
     # initialization
     randomize <- NULL
     plan <- list()
     
     for (i in 1:trial) {
          if (case == 1) {
               # the following code was revised by AAGulles c/o VIBartolome 27May 2014
               # the argument BlockIn2D was removed 
               # and replaced by the arguments RowsInBlockSequence and ColumnsInBlockSequence
               result <- try(temp <- DiGGer(NumberOfTreatments = length(tempComb[[1]]),
                                            RowsInDesign = blksize * r,
                                            ColumnsInDesign = numBlk,
                                            RowsInReplicate = blksize,
                                            ColumnsInReplicate = numBlk,
                                            # BlockIn2D =c(blksize * r, 1),
                                            RowsInBlockSequence = c(blksize * r, blksize),
                                            ColumnsInBlockSequence = c(numBlk, 1),
                                            TreatmentName = tempComb[[1]]),
                             silent = TRUE)
          } else {
               # the following code was revised by AAGulles c/o VIBartolome 27May 2014
               # the argument BlockIn2D was removed 
               # and replaced by the arguments RowsInBlockSequence and ColumnsInBlockSequence
               result <- try(temp <- DiGGer(NumberOfTreatments = length(tempComb[[1]]),
                                            RowsInDesign = numBlk,
                                            ColumnsInDesign = blksize * r,
                                            RowsInReplicate = numBlk,
                                            ColumnsInReplicate = blksize,
                                            #BlockIn2D =c(1, blksize * r),
                                            RowsInBlockSequence = c(numBlk, 1),
                                            ColumnsInBlockSequence = c(blksize*r, blksize),
                                            TreatmentName = tempComb[[1]]),
                             silent = TRUE)
          }
          
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          capture.output(run(temp))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          #plot(temp)
          result_mat <- getDesign(temp)   # get the layout
          fbook <- print(temp, option = "list") #  save the fieldbook of the final design to an R object
          plan[[i]] <- matrix(fbook$ID, nrow(result_mat), ncol(result_mat)) # get the layout with the specified treatment level
          if (i == 1) {
               if (case == 1) {
                    blkNum <- matrix(1:numBlk, nrow(plan[[i]]), ncol(plan[[i]]), byrow = TRUE)
                    repNum <- matrix(as.numeric(paste(fbook$REP,paste(c(rep(0, max(nchar(1:length(tempComb[[1]]))))), collapse = ""), sep = "")), nrow(plan[[i]]), ncol(plan[[i]]))
                    plotNumTemp <- matrix(1:length(tempComb[[1]]), blksize, numBlk, byrow = TRUE)
                    if (serpentine) { for (k in seq(2, blksize, by = 2)) { plotNumTemp[k,] <- rev(plotNumTemp[k,]) }}
                    plotNum <- NULL
                    for (j in 1:r) { plotNum <- rbind(plotNum, plotNumTemp) }
                    plotNum <- repNum + plotNum
               } else {
                    blkNum <- matrix(1:numBlk, nrow(plan[[i]]), ncol(plan[[i]]), byrow = FALSE)
                    repNum <- matrix(as.numeric(paste(fbook$REP,paste(c(rep(0, max(nchar(1:length(tempComb[[1]]))))), collapse = ""), sep = "")), nrow(plan[[i]]), ncol(plan[[i]]))
                    plotNumTemp <- matrix(1:length(tempComb[[1]]), numBlk, blksize, byrow = TRUE)
                    if (serpentine) { for (k in seq(2, numBlk, by = 2)) { plotNumTemp[k,] <- rev(plotNumTemp[k,]) }}
                    plotNum <- NULL
                    for (j in 1:r) { plotNum <- cbind(plotNum, plotNumTemp) }
                    plotNum <- repNum + plotNum
               }
               
          }
          
          tempFieldOrder <- merge(as.data.frame.table(blkNum), as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("Block","PlotNum"))
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"]) 
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          colnames(tempFieldOrder) <- gsub("Freq", "", colnames(tempFieldOrder))
          randomize <- rbind(randomize, cbind(Trial = 1, merge(print(temp, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
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
     cat("\t","Latinized Alpha Lattice Design","\n\n",sep = "") 
     cat(toupper("Design Parameters:"),"\n",sep = "")
     cat("\t","Number of Trials = ", trial, "\n",sep = "")
     cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
     cat("\t","Number of Replicates = ", r, "\n",sep = "")
     cat("\t","Number of Plots per Block = ", blksize, "\n",sep = "")
     cat("\t","Number of Blocks per Replicate = ", length(tempComb[[1]])/blksize, "\n\n",sep = "")
     
     cat("\t","Number of Field Rows = ", numFieldRow, "\n\n",sep = "")
     #cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")
     
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
     
     
} # end of designLatinizedAlpha.default

