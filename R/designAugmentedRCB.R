# -------------------------------------------------------------------------------------
# designAugmentedRCB: Generate randomization for augmented design in RCBD.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 11.27.2014
# Notes: Revised the old version. 
# Features: Equal number of rows per block
#           Block layout can be determine by the user
# -------------------------------------------------------------------------------------

designAugmentedRCB <- function(numCheck, numNew, trmtName = NULL, r = 2, trial = 1, 
                               rowPerBlk, numFieldRow, serpentine = FALSE, 
                               trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL, 
                               file = NULL) UseMethod("designAugmentedRCB")

designAugmentedRCB.default <- function(numCheck, numNew, trmtName = NULL, r = 2, trial = 1, 
                                       rowPerBlk, numFieldRow, serpentine = FALSE, 
                                       trmtLabel = NULL, checkTrmt = NULL, newTrmt = NULL, 
                                       file = NULL) {
     
     # check user input
     # -- check the number of replicates
     if (r < 2) { stop("The number of replicates should be greater than or equal to 2.") }
     
     # -- determine if the number of elements in checkTrmt/newTrmt is equal to numCheck/numNew
     withCheckTrmtEntry <- FALSE
     withNewTrmtEntry <- FALSE
     if (is.null(checkTrmt)) { checkTrmt <- paste("check",1:numCheck, sep = "")
     } else {
          withCheckTrmtEntry <- TRUE
          if (length(checkTrmt) != numCheck) { stop("Error: The number of elements of the arg 'checkTrmt' is not equal to numCheck.") } 
          #checkTrmt <- as.character(sample(as.character(checkTrmt), numCheck))
     }
     
     if (is.null(newTrmt)) { newTrmt <- paste("new",1:numNew, sep = "")
     } else {
          withNewTrmtEntry <- TRUE
          if (length(newTrmt) != numNew) { stop("Error: The number of elements of the arg 'newTrmt' is not equal to numNew.") }          
          #newTrmt <- as.character(sample(as.character(newTrmt), numNew))
     }
     
     # --- randomize the checkTrmt and newTrmt
     #checkTrmt <- sample(checkTrmt, numCheck, replace = FALSE)
     #newTrmt <- sample(newTrmt, numNew, replace = FALSE)
     
     numTrmt <- numCheck + numNew            # -- determine the total number of treatment in the experiment
     numNewPerBlk <- numNew/r                # -- determine the number of newTrmt per blk
     numPlots <- (numCheck * r) + numNew     # -- determine the total number of experimental units 
     errorDF <- (r - 1) * (length(checkTrmt) - 1)

     #if (numPlots %% numFieldRow != 0)
     numFieldCol <- numPlots/numFieldRow     # -- determine the number of columns in the design
     
     #if (numPlots %% r != 0)
     numPlotPerBlk <- numPlots/r            # -- assume equal number of plots per block
     colPerBlk <- numPlotPerBlk/rowPerBlk   # -- determine the number of columns per block
     
     # -- determine the orientation of blocks
     
     numBlkRow <- numFieldRow/rowPerBlk
     numBlkCol <- numFieldCol/colPerBlk
     if (numBlkRow * numBlkCol != r) { stop("?Number of field rows not divisible by the number of rows per block.") }
     
     if (!is.null(file)) { 
          tempFile <- strsplit(file, split = "\\.")[[1]]
          tempExt <- tolower(tempFile[length(tempFile)])
          if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
          newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
     }
     
     # generate randomization and fieldbook
     capture.output(result <- designRCBD(generate = list(EntryNo = 1:numPlotPerBlk), r, trial, numFieldRow, rowPerBlk, serpentine, display = FALSE))
     
     fbook <- result$fieldbook
     if (is.factor(fbook[,"EntryNo"])) { fbook[,"EntryNo"] <- as.numeric(fbook[,"EntryNo"]) }
     #fbook[,"ID"] <- as.character(fbook[,"EntryNo"])
     fbook[,"ID"] <- NA
     
     #-- recoding the entry number for augmented design in RCB
     checkEntryNoOld <- (numPlotPerBlk-numCheck+1):numPlotPerBlk
     checkEntryNoNew <- sample((numTrmt-numCheck+1):numTrmt, numCheck)
     for (i in (1:length(checkEntryNoOld))) { 
          fbook[fbook[,"EntryNo"] == checkEntryNoOld[i], "ID"] <- checkTrmt[i]
          fbook[fbook[,"EntryNo"] == checkEntryNoOld[i], "EntryNo"] <- checkEntryNoNew[i] 
     }
     newEntryNoOld <- 1:numNewPerBlk
     if (withNewTrmtEntry) { newEntryNoNew <- sample(1:numNew, numNew)        
     } else { newEntryNoNew <- 1:numNew }
     #newEntryNoNew <- sample(1:numNew, numNew)
     if (trial > 1) { for (i in (2:trial)) { 
          newEntryNoNew <- c(newEntryNoNew, newEntryNoNew) 
          newTrmt <- c(newTrmt, newTrmt)
     }}
     fbook <- fbook[order(fbook$Trial, fbook$EntryNo),]
     fbook[fbook[,"EntryNo"] %in% newEntryNoOld, "ID"] <- newTrmt
     fbook[fbook[,"EntryNo"] %in% newEntryNoOld, "EntryNo"] <- newEntryNoNew
     fbook[,"EntryNo"] <- as.factor(fbook[,"EntryNo"]) 
     fbook[,"ID"] <- as.factor(fbook[,"ID"])
     
     #-- construct the new layout for the augmented design in RCB
     fbook <- fbook[order(fbook$Trial, fbook$FieldRow, fbook$FieldCol),]
     result$layout1 <- list()
     for (i in (1:trial)) {
          result$layout[[i]] <- matrix(fbook[fbook[,"Trial"] == i,"EntryNo"],
                                       nrow = nrow(result$plotNum), 
                                       ncol = ncol(result$plotNum), 
                                       byrow = TRUE,
                                       dimnames = dimnames(result$plotNum))

          result$layout1[[i]] <- matrix(fbook[fbook[,"Trial"] == i,"ID"],
                                       nrow = nrow(result$plotNum), 
                                       ncol = ncol(result$plotNum), 
                                       byrow = TRUE,
                                       dimnames = dimnames(result$plotNum))
          
          
     }
     fbook <- fbook[,c("Trial", "Rep", "ID", "EntryNo", "PlotNum", "FieldRow", "FieldCol")]
     
     # rename the fbook if user supply a trmtName and trmtLabel
     if (!is.null(trmtName))  { names(fbook)[match("EntryNo", names(fbook))] <- trmtName  }
     if (!is.null(trmtLabel)) { names(fbook)[match("ID", names(fbook))] <- trmtLabel  }
     
     result$fieldbook <- fbook[order(fbook$Trial, fbook$PlotNum),]
     
     # printing of the layout
     for (i in (1:trial)) {
          if (!is.null(file)) { png(filename = paste(newFile, "Layout_Trial",i,".png",sep = ""))  }
          des.plot(result$layout[[i]], seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
                   bwd = 4, bcol = 4, cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
          des.plot(result$layout[[i]], seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
                   bdef = cbind(rowPerBlk, colPerBlk), bwd = 4)
          #des.plot(result$layout[[i]], seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3)
          #des.plot(result$layout[[i]], bdef = cbind(rowPerBlk, colPerBlk), col = 7, bwd = 4, new = FALSE, chtdiv = 3)
          if (!is.null(file)) { dev.off() }  
     }
     
     cat(toupper("Design Properties:"), "\n", sep = "")
     cat("\t", "Augmented Randomized Complete Block Design (Augmented RCBD)", "\n\n", sep = "")
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
     if (errorDF < 12) { cat("WARNING: Too few error df.","\n\n") }
     
     # saving the fieldbook to a file
     if (!is.null(file)) {
          if (tempExt == "csv") { write.csv(result$fieldbook, file = paste(newFile, tempExt, sep = "."), row.names = FALSE)
          } else { save(result$fieldbook, file = paste(newFile, tempExt, sep = ".")) }
     } else {
          cat("Results of Randomization: \n")
          printDataFrame(result$fieldbook)
     }
   
     return(result)
}
