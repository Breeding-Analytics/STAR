# -------------------------------------------------------------------------------------
# designAugmentedRowColumn
# Description: Generate randomization and layout for augmented Row-Column design.
# Script Created by: Violeta I Bartolome for International Rice Research Institute
# Function Created by: Alaine A. Gulles 10.15.2014 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.08.2014
# Note: Uses the DiGGer function from package DiGGer
# -------------------------------------------------------------------------------------
# Arguments:
# numCheck
# -------------------------------------------------------------------------------------

designAugmentedRowColumn <- function(numCheck, numNew, trmtName = NULL, r = 2, trial = 1, 
                                     rowblkPerRep, rowPerRowblk, numFieldRow, serpentine = FALSE, 
                                     trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL,
                                     file = NULL) UseMethod("designAugmentedRowColumn")

designAugmentedRowColumn.default <- function(numCheck, numNew, trmtName = NULL, r = 2, trial = 1, 
                                     rowblkPerRep, rowPerRowblk, numFieldRow, serpentine = FALSE, 
                                     trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL,
                                     file = NULL) {
     # check user input
     # -- check the number of replicates
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.") }
       
     # -- determine if the number of elements in checkTrmt/newTrmt is equal to numCheck/numNew
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
     checkTrmt <- sample(checkTrmt, numCheck, replace = FALSE)
     newTrmt <- sample(newTrmt, numNew, replace = FALSE)
     
     # --- determine the total number of experimental units 
     numTrmt <- numCheck + numNew 
     numPlots <- (numCheck * r) + numNew
     
     if (numPlots > 1500) { stop("The maximum number of experimental units that can be generated is 1500.") }
     if (numPlots%%numFieldRow != 0) { stop("Total number of experimental units should be divisible by the number of field rows.") }
     if (numFieldRow%%(rowblkPerRep*rowPerRowblk) != 0) { stop("The number of fieldrows should be divisible by the product of the number of row block per rep and the number of rows per row block.") } 
     numFieldCol <- numPlots/numFieldRow
     
     if (numPlots%%r != 0) { stop("The total number of experimental units should be divisible by the number of replicates.") }
     numPlotPerRep <- numPlots/r
     
     rowPerRep <- rowblkPerRep * rowPerRowblk # total number of row per rep
     
     colblkPerRep <- numCheck/rowblkPerRep
     
     #colPerRep <- colblkPerRep * colPerColblk
     colPerRep <- numPlotPerRep/rowPerRep
     
     #colPerColblk <- numFieldCol/colblkPerRep
     colPerColblk <-  colPerRep/colblkPerRep
     
     
     
     
     numRepRow <- numFieldRow/rowPerRep 
     numRepCol <- numFieldCol/colPerRep                     # determine the number of blocks along the length of each replicate
     
     randomize <- NULL
     plan <- list()
     plan1 <- list()
     plotNum <- NULL
     blkNum <- NULL
     repNum <- NULL

     if (!is.null(file)) { 
          tempFile <- strsplit(file, split = "\\.")[[1]]
          tempExt <- tolower(tempFile[length(tempFile)])
          if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
          newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
     }
     
     for (i in (1:trial)) {
          capture.output(result <- try(augRC <- des.prep00(nt = numTrmt,
                                            nrd = numFieldRow,
                                            ncd = numFieldCol,
                                            trep = rep(c(1, r), c(numNew, numCheck)),
                                            tgrp = rep(c(1, 2), c(numNew, numCheck)),
                                            ribs = c(rowPerRep, rowPerRowblk),
                                            cibs = c(colPerRep, colPerColblk),
                                            tnam = c(newTrmt, checkTrmt),
                                            tnum = 1:numTrmt),
                        silent = TRUE))
          
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          capture.output(run(augRC))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          augRC_mat <- getDesign(augRC)
          if (!is.null(file)) { png(filename = paste(newFile, "Layout_Trial",i,".png",sep = ""))  }
          #des.plot(augRC_mat, seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
          #         bdef = cbind(rowblkPerRep, colblkPerRep), bwd = 4, bcol = 4,
          #         cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
          #des.plot(augRC_mat, seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
          #         bdef = cbind(rowPerRep, colPerRep), bwd = 4)
          des.plot(augRC_mat, seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3, bwd = 4, bcol = 4,
                   cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
          des.plot(augRC_mat, seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
                   bdef = cbind(rowPerRowblk, colPerColblk), bwd = 3)
          if (!is.null(file)) { dev.off() }  
          
          plan[[i]] <- matrix(print(augRC, option = "list")$ENTRY, nrow(augRC_mat), ncol(augRC_mat))
          plan1[[i]] <- matrix(print(augRC, option = "list")$ID, nrow(augRC_mat), ncol(augRC_mat))
          
          if (i == 1) {
               
               rowblkNumTemp <- matrix(0, rowPerRep, colPerRep)
               colblkNumTemp <- matrix(0, rowPerRep, colPerRep)
               
               tempPlotNumRep <- matrix(1:numPlotPerRep, rowPerRep, colPerRep, byrow = TRUE)
               if (serpentine) { for (k in seq(2, rowPerRep, by = 2)) { tempPlotNumRep[k,] <- rev(tempPlotNumRep[k,]) }}
               
               for (j in 1:rowblkPerRep) {
                    rowIndexLL <- (j * rowPerRowblk) - rowPerRowblk + 1
                    rowIndexUL <- rowIndexLL + rowPerRowblk - 1
                    for (k in 1:colblkPerRep) {
                         colIndexLL <- (k * colPerColblk) - colPerColblk + 1
                         colIndexUL <- colIndexLL + colPerColblk - 1
                         rowblkNumTemp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <-  j
                         colblkNumTemp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <-  k
                    } ## end stmt -- for (k in 1:numBlkCol)
               } ## end stmt -- for (j in 1:numBlkRow)
               
               repNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               plotNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               rowblkNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               colblkNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               
               repcode <- 1
               for (j in 1:numRepRow) {
                    rowIndexLL <- (j * rowPerRep) - rowPerRep + 1
                    rowIndexUL <- rowIndexLL + rowPerRep - 1
                    for (k in 1:numRepCol) {
                         colIndexLL <- (k * colPerRep) - colPerRep + 1
                         colIndexUL <- colIndexLL + colPerRep - 1
                         repNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- repcode
                         rowblkNum [rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- rowblkNumTemp
                         colblkNum [rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- colblkNumTemp
                         plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- as.numeric(paste(repcode,paste(c(rep(0, nchar(numPlotPerRep))), collapse = ""), sep = "")) + tempPlotNumRep
                         repcode <- repcode + 1 
                    } ## end stmt -- for (k in 1:numRepCol)
               } ## end stmt -- for (j in 1:numRepRow)
          } ## end stmt -- if (i == 1)
          
          tempFieldOrder <- merge(as.data.frame.table(repNum), as.data.frame.table(rowblkNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("REP", "ROWBLK")) 
          tempFieldOrder <- merge(tempFieldOrder, as.data.frame.table(colblkNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
          tempFieldOrder <- merge(tempFieldOrder, as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("COLBLK", "PLOTNO"))
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"])
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          randomize <- rbind(randomize, cbind(TRIAL = i, merge(print(augRC, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
          dimnames(plan[[i]]) <- list(paste("FIELDROW", 1:nrow(plan[[i]]), sep = ""), paste("FIELDCOL", 1:ncol(plan[[i]]), sep = ""))
          dimnames(plan1[[i]]) <- dimnames(plan[[i]])
     }  ## end stmt -- for (i in 1:trial)
     
     dimnames(plotNum) <- dimnames(plan[[1]])
     dimnames(rowblkNum) <- dimnames(plan[[1]])
     dimnames(colblkNum) <- dimnames(plan[[1]])
     dimnames(repNum) <- dimnames(plan[[1]])
     names(plan) <- paste("Trial", 1:trial, sep = "")
     names(plan1) <- names(plan)
     randomize <- randomize[,c("TRIAL", "FreqREP", "FreqROWBLK", "FreqCOLBLK", "ID", "ENTRY", "FreqPLOTNO", "ROW", "RANGE")]
     names(randomize) <- c("Trial", "Rep", "RowBlk", "ColBlk", "ID", "ENTRY", "PlotNum", "FieldRow", "FieldCol")
     #colnames(randomize) <- gsub("Freq","", colnames(randomize))
     #colnames(randomize)[match("ROW", names(randomize))] <- "FieldRow"
     #colnames(randomize)[match("RANGE", names(randomize))] <- "FieldCol"
     if (!is.null(trmtName)) { colnames(randomize)[match("ENTRY", names(randomize))] <- trmtName
     } else { colnames(randomize)[match("ENTRY", names(randomize))] <- "EntryNo" }
     if (!is.null(trmtLabel)) { colnames(randomize)[match("ID", names(randomize))] <- trmtLabel }
     
     # re-arrange the fieldbook
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum), ]
     rownames(randomize) <- 1:nrow(randomize)
     
     # print of design parameters:
     cat(toupper("Design Properties:"), "\n", sep = "")
     cat("\t", "Incomplete Block Design", "\n", sep = "")
     cat("\t", "Augmented Row-Column Design", "\n\n", sep = "")
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
          cat(", ..., ", newTrmt[numNew], "\n\n", sep = "") 
     }
     cat("\t", "Number of Field Rows = ", numFieldRow,"\n", sep = "")
     cat("\t", "Number of Field Columns = ", numFieldCol, "\n\n", sep = "")
     
     if (!is.null(file)) {
          if (tempExt == "csv") { write.csv(randomize, file = paste(newFile, tempExt, sep = "."), row.names = FALSE)
          } else { save(randomize, file = paste(newFile, tempExt, sep = "."), row.names = FALSE) }
     } else {
          cat("Results of Randomization:\n")
          printDataFrame(randomize)
     }
     
     if (is.null(trmtLabel)) { return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum, repNum = repNum)))
     } else { return(invisible(list(fieldbook = randomize, layout = plan, layoutID = plan1, plotNum = plotNum, repNum = repNum))) }
} ## end stmt -- designAugmentedRowColumn.default