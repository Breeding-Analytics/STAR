# -------------------------------------------------------------------------------------
# designAugmentedAlpha
# Description: Generate randomization and layout for augmented design in Alpha Lattice.
# Script Created by: Violeta I Bartolome for International Rice Research Institute
# Function Created by: Alaine A. Gulles 05.29.2014 for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.29.2014
# Note: Uses the DiGGer function from package DiGGer
# -------------------------------------------------------------------------------------
# Arguments:
# numCheck
# -------------------------------------------------------------------------------------

designAugmentedAlpha <- function(numCheck, numNew, trmtName = NULL, blksize, r = 2, trial = 1, 
                                 rowPerBlk, rowPerRep, numFieldRow, 
                                 serpentine = FALSE, trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL,
                                 file = NULL) UseMethod("designAugmentedAlpha")

designAugmentedAlpha.default <- function(numCheck, numNew, trmtName = NULL, blksize, r = 2, trial = 1, 
                                         rowPerBlk, rowPerRep, numFieldRow, 
                                         serpentine = FALSE, trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL,
                                         file = NULL) {

     # Check input of the user:
     # --- check of the number of replicate is at least 2
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.") }
     if (r > 4) { stop("Tne maximum number of replicates that can be generated for this design is 4.") }
     
     # --- determine if the number of elements in checkTrmt/newTrmt is equal to numCheck/numNew
     if (is.null(checkTrmt)) { checkTrmt <- paste("check",1:numCheck, sep = "")
     } else {
          if (length(checkTrmt) != numCheck) { stop("Error: The number of elements of the arg 'checkTrmt' is not equal to numCheck.") }          
     }
     if (is.null(newTrmt)) { newTrmt <- paste("new",1:numNew, sep = "")
     } else {
          if (length(newTrmt) != numNew) { stop("Error: The number of elements of the arg 'newTrmt' is not equal to numNew.") }          
     }
     checkTrmt <- as.character(checkTrmt)
     newTrmt <- as.character(newTrmt)
     
     # --- randomize the checkTrmt and newTrmt
     #checkTrmt <- sample(checkTrmt, numCheck, replace = FALSE)
     #newTrmt <- sample(newTrmt, numNew, replace = FALSE)
     
     # determine the total number of experimental units    
     numTrmt <- numCheck + numNew 
     numPlots <- (numCheck * r) + numNew
     
     if (numPlots > 1500) { stop("The maximum number of experimental units that can be generated is 1500.") }
     if (numPlots %% numFieldRow != 0) { stop("Total number of experimental units should be divisible by the number of field rows.") }
     numFieldCol <- numPlots/numFieldRow
     
     if (numPlots%%r != 0) { stop("The total number of experimental units should be divisible by the number of replicates.") }
     numPlotPerRep <- numPlots/r
     
     if (numPlotPerRep%%rowPerRep != 0) { stop("The total number of experimental units per replicated should be divisible by the number of rows per replicates.") }
     colPerRep <- numPlotPerRep/rowPerRep
     
     if (blksize%%rowPerBlk != 0) { stop("The total number of plots per block should be divisible by the number of rows per block.") }
     colPerBlk <- blksize/rowPerBlk
     
     numRepRow <- numFieldRow/rowPerRep 
     numRepCol <- numFieldCol/colPerRep                     # determine the number of blocks along the length of each replicate
     numBlkRow <- rowPerRep/rowPerBlk                       # determine the number of blocks along the length of each replicate
     numBlkCol <- colPerRep/colPerBlk 
     
     randomize <- NULL
     plan <- list()
     plan1 <- list()
     plotNum <- NULL
     blkNum <- NULL
     repNum <- NULL
          
     if (rowPerBlk == 1) serpentine <- FALSE
     
     if (!is.null(file)) { 
          tempFile <- strsplit(file, split = "\\.")[[1]]
          tempExt <- tolower(tempFile[length(tempFile)])
          if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
          newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
     }
     
     for (i in 1:trial) {
          capture.output(result <- try(augAlpha <- des.prep00(nt = numTrmt,
                                               nrd = numFieldRow,
                                               ncd = numFieldCol,
                                               trep = rep(c(1, r), c(numNew, numCheck)),
                                               tgrp = rep(c(1, 2), c(numNew, numCheck)),
                                               ribs = c(rowPerRep, rowPerBlk),
                                               cibs = c(colPerRep, colPerBlk),
                                               tnam = c(newTrmt, checkTrmt),
                                               tnum = 1:numTrmt),
                        silent = TRUE))
          
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          capture.output(run(augAlpha))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          augAlpha_mat <- getDesign(augAlpha)
          if (!is.null(file)) { png(filename = paste(newFile, "Layout_Trial",i,".png",sep = ""))  }
          #des.plot(augAlpha_mat, seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
          #         bdef = cbind(rowPerBlk, colPerBlk), bwd = 4, bcol = 4,
          #         cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
          #des.plot(augAlpha_mat, seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
          #         bdef = cbind(rowPerRep, colPerRep), bwd = 4)
          
          
          des.plot(augAlpha_mat, seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
                   cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
          des.plot(augAlpha_mat, seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
                   bdef = cbind(rowPerBlk, colPerBlk), bwd = 4, bcol = 4)
          des.plot(augAlpha_mat, NULL, new = FALSE, label = TRUE, chtdiv = 3,
                   bdef = cbind(rowPerRep, colPerRep), bwd = 4)
          
          if (!is.null(file)) { dev.off() }           
          
          plan[[i]] <- matrix(print(augAlpha, option = "list")$ENTRY, nrow(augAlpha_mat), ncol(augAlpha_mat))
          plan1[[i]] <- matrix(print(augAlpha, option = "list")$ID, nrow(augAlpha_mat), ncol(augAlpha_mat))
          
          if (i == 1) {
               
               blkNumTemp <- matrix(0, rowPerRep, colPerRep)
               blkPlotNumTemp <- matrix(0, rowPerRep, colPerRep)
               
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
                    } ## end stmt -- for (k in 1:numBlkCol)
               } ## end stmt -- for (j in 1:numBlkRow)
               
               
               blkNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               repNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               plotNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               
               repcode <- 1
               for (j in 1:numRepRow) {
                    for (k in 1:numRepCol) {
                         rowIndexLL <- (j * rowPerRep) - rowPerRep + 1
                         rowIndexUL <- rowIndexLL + rowPerRep - 1
                         colIndexLL <- (k * colPerRep) - colPerRep + 1
                         colIndexUL <- colIndexLL + colPerRep - 1
                         repNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- repcode
                         plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- as.numeric(paste(repcode,paste(c(rep(0, nchar(numTrmt))), collapse = ""), sep = "")) + blkPlotNumTemp
                         blkNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- blkNumTemp
                         repcode <- repcode + 1 
                    } ## end stmt -- for (k in 1:numRepCol)
               } ## end stmt -- for (j in 1:numRepRow)
          } ## end stmt -- if (i == 1)
          
          tempFieldOrder <- merge(as.data.frame.table(repNum), as.data.frame.table(blkNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("REP", "BLK")) 
          tempFieldOrder <- merge(tempFieldOrder, as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"])
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          randomize <- rbind(randomize, cbind(Trial = i, merge(print(augAlpha, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
          dimnames(plan[[i]]) <- list(paste("FIELDROW", 1:nrow(plan[[i]]), sep = ""), paste("FIELDCOL", 1:ncol(plan[[i]]), sep = ""))
          dimnames(plan1[[i]]) <- dimnames(plan[[i]])
     }  ## end stmt -- for (i in 1:trial)
     
     dimnames(plotNum) <- dimnames(plan[[1]])
     dimnames(blkNum) <- dimnames(plan[[1]])
     dimnames(repNum) <- dimnames(plan[[1]])
     names(plan) <- paste("Trial", 1:trial, sep = "")
     names(plan1) <- names(plan)
     randomize <- randomize[,c("Trial", "FreqREP", "FreqBLK", "ID", "ENTRY", "Freq", "ROW", "RANGE")]
     colnames(randomize)[match("Freq", names(randomize))] <- "PlotNum"
     #colnames(randomize) <- gsub("Freq","", colnames(randomize))
     colnames(randomize)[match("FreqREP", names(randomize))] <- "Rep"
     colnames(randomize)[match("FreqBLK", names(randomize))] <- "Block"
     colnames(randomize)[match("ROW", names(randomize))] <- "FieldRow"
     colnames(randomize)[match("RANGE", names(randomize))] <- "FieldCol"
     if (!is.null(trmtName)) { colnames(randomize)[match("ENTRY", names(randomize))] <- trmtName
     } else { colnames(randomize)[match("ENTRY", names(randomize))] <- "EntryNo" }
     if (!is.null(trmtLabel)) { colnames(randomize)[match("ID", names(randomize))] <- trmtLabel }
     
     # re-arrange the fieldbook
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum), ]
     rownames(randomize) <- 1:nrow(randomize)
     
     # print of design parameters:
     cat(toupper("Design Properties:"), "\n", sep = "")
     cat("\t", "Incomplete Block Design", "\n", sep = "")
     cat("\t", "Augmented Alpha Lattice Design", "\n\n", sep = "")
     cat(toupper("Design Parameters"), "\n", sep = "")
     cat("\t", "Number of Trials = ", trial, "\n", sep = "")
     cat("\t", "Number of Replicated Treatments = ", numCheck, "\n", sep = "")
     cat("\t", "Levels of Replicated Treatments = ", sep = "")
     if (numCheck <= 5) { cat(paste(checkTrmt, collapse = ", ", sep = ""), "\n", sep = "") 
     } else {
          cat(paste(checkTrmt[1:3], collapse = ", ", sep = ""), sep = "") 
          cat(", ..., ", checkTrmt[numCheck], "\n", sep = "") 
     }
     cat("\t", "Number of Replicates = ", r, "\n", sep = "")
     cat("\t", "Number of UnReplicated Treatments = ", numNew, "\n", sep = "")
     cat("\t", "Levels of UnReplicated Treatments = ", sep = "")
     if (numNew <= 5) { cat(paste(newTrmt, collapse = ", ", sep = ""), "\n", sep = "") 
     } else {
          cat(paste(newTrmt[1:3], collapse = ", ", sep = ""), sep = "") 
          cat(", ..., ", newTrmt[numNew], "\n", sep = "") 
     }
     cat("\t", "Number of Plots per Block = ", blksize,"\n", sep = "")
     cat("\t", "Number of Blocks per Replicate = ", numPlotPerRep/blksize,"\n\n", sep = "")
     cat("\t", "Number of Field Rows = ", numFieldRow,"\n", sep = "")
     cat("\t", "Number of Field Columns = ", numFieldCol, "\n\n", sep = "")
     
     if (!is.null(file)) {
          if (tempExt == "csv") { write.csv(randomize, file = paste(newFile, tempExt, sep = "."), row.names = FALSE)
          } else { save(randomize, file = paste(newFile, tempExt, sep = "."), row.names = FALSE) }
     } else {
          cat("Results of Randomization:\n")
          printDataFrame(randomize)
     }
     
     if (is.null(trmtLabel)) { return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum, blkNum = blkNum, repNum = repNum)))
     } else { return(invisible(list(fieldbook = randomize, layout = plan, layoutID = plan1, plotNum = plotNum, blkNum = blkNum, repNum = repNum))) }
} ## end stmt -- designAugmentedAlpha

