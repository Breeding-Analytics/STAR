# ------------------------------------------------------------------------------------------
# designAugmentedLSD: Generate randomization for augmented design in Latin Square Design.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.10.2013
# ------------------------------------------------------------------------------------------

designAugmentedLSD <- function(numCheck, numNew, trmtName = NULL, trial = 1, 
                               numFieldRow, serpentine = FALSE, trmtLabel = NULL, 
                               checkTrmt = NULL, newTrmt = NULL, file = NULL) UseMethod("designAugmentedLSD")

designAugmentedLSD.default <- function(numCheck, numNew, trmtName = NULL, trial = 1, 
                                       numFieldRow, serpentine = FALSE, trmtLabel = NULL, 
                                       checkTrmt = NULL, newTrmt = NULL, file = NULL) {

     # -- determine if the number of elements in checkTrmt/newTrmt is equal to numCheck/numNew
     withCheckTrmtEntry <- FALSE
     withNewTrmtEntry <- FALSE
     if (is.null(checkTrmt)) { checkTrmt <- paste("check",1:numCheck, sep = "")
     } else {
          withCheckTrmtEntry <- TRUE
          if (length(checkTrmt) != numCheck) { stop("Error: The number of elements of the arg 'checkTrmt' is not equal to numCheck.") } 
     }
     
     if (is.null(newTrmt)) { newTrmt <- paste("new",1:numNew, sep = "")
     } else {
          withNewTrmtEntry <- TRUE
          if (length(newTrmt) != numNew) { stop("Error: The number of elements of the arg 'newTrmt' is not equal to numNew.") }          
     }

     # CHECK
     if (length(checkTrmt) > numFieldRow) { stop("Too few rows for blocking.") }
     totalNumPlots <- length(checkTrmt)**2 + length(newTrmt)
     if (totalNumPlots%%numFieldRow != 0) { stop("Number of field rows must divide number of plots.") } 
     numFieldCol <- totalNumPlots/numFieldRow
     if (length(checkTrmt) > numFieldCol) { stop("Too few columns for blocking.") }
     errorDf <- (length(checkTrmt) - 1)*(length(checkTrmt) - 2)
     
     embedColInBasicCol <- floor(numFieldCol/length(checkTrmt))
     basicColWAddlCol <- numFieldCol - (length(checkTrmt) * embedColInBasicCol)
     
     plotsInBasicRow <- totalNumPlots/length(checkTrmt) # -- number of plots within basic row of LSD
     embedRowInBasicRow <- floor(numFieldRow/length(checkTrmt)) # -- number of embedded row within the basic row of LSD 
     basicRowWAddlRow <- numFieldRow - (length(checkTrmt) *  embedRowInBasicRow)
     
     plotsInFieldRow <- totalNumPlots/numFieldRow
     plotsInFieldCol <- totalNumPlots/numFieldCol
     
     EntryCodeNew <- 1:numNew
     EntryCodeCheck <- (numNew+1):(numCheck + numNew)
          
     randomize <- NULL
     plan <- list()
     plotNum <- list()
     capture.output(result <- designLSD(list(tempTrmt = EntryCodeCheck), trial, serpentine = FALSE, display = FALSE))
     tempDesign <- result$fieldbook
     for (i in (1:trial)) {
          tempSample <- tempDesign[tempDesign[,"Trial"] == 1,]
          tempPlan <- result$layout[[i]]  
          tempNewDesign <- NULL
          tempColDesign <- NULL
          if (embedRowInBasicRow == 1 && basicRowWAddlRow == 0) { tempNewDesign <- tempPlan
          } else {
               for (j in (1:length(checkTrmt))) {
                    if (j <= basicRowWAddlRow) { tempRow <- matrix(0, nrow = (embedRowInBasicRow + 1), ncol = length(checkTrmt))
                    } else { tempRow <- matrix(0, nrow = embedRowInBasicRow, ncol = length(checkTrmt)) }
                    for (k in (1:length(checkTrmt))) { 
                         if (j <= basicRowWAddlRow) { tempRow[sample(1:(embedRowInBasicRow + 1),1),k] <- tempPlan[j,k] 
                         } else { tempRow[sample(1:embedRowInBasicRow,1),k] <- tempPlan[j,k] }
                    }
                    tempNewDesign <- rbind(tempNewDesign, tempRow)
               }                
          }
          
          if (embedColInBasicCol == 1 && basicColWAddlCol == 0) { tempColDesign <- tempNewDesign
          } else {
               for (k in (1:length(checkTrmt))) {
                    if (k <= basicColWAddlCol) { tempCol <- matrix(0, nrow = nrow(tempNewDesign), ncol = (embedColInBasicCol + 1))
                    } else { tempCol <- matrix(0, nrow = nrow(tempNewDesign), ncol = embedColInBasicCol) }
                    for (j in (1:nrow(tempNewDesign))) { 
                         if (k <= basicColWAddlCol) { tempCol[j,sample(1:(embedColInBasicCol + 1),1)] <- tempNewDesign[j,k] 
                         } else { tempCol[j,sample(1:embedColInBasicCol,1)] <- tempNewDesign[j,k]  }
                    }
                    tempColDesign <- cbind(tempColDesign, tempCol)
               }               
          }
          #tempColDesign[tempColDesign == "0"] <- sample(newTrmt, length(newTrmt))
          tempColDesign[tempColDesign == "0"] <- sample(1:numNew, length(newTrmt))
          tempNewDesign <- as.data.frame.table(tempColDesign)
          if (basicRowWAddlRow == 0) { tempNewDesign$Row <- rep(1:length(checkTrmt), each = embedRowInBasicRow)
          } else {  
               tempNewDesign$Row <- c(rep(1:basicRowWAddlRow, each = (embedRowInBasicRow + 1)), 
                                      rep((basicRowWAddlRow+1):length(checkTrmt), each = embedRowInBasicRow))
          }
          if (basicColWAddlCol == 0) { tempNewDesign$Column <- rep(1:length(checkTrmt), each = nrow(tempColDesign)*embedColInBasicCol)
          } else {
               tempNewDesign$Column <- c(rep(1:basicColWAddlCol, each = ((embedColInBasicCol + 1)*nrow(tempColDesign))), 
                                         rep((basicColWAddlCol+1):length(checkTrmt), each = (embedColInBasicCol*nrow(tempColDesign))))
          }
          tempNewDesign$Trial <- i
          plan[[i]] <- tempColDesign
          plotNum[[i]] <- matrix(1:totalNumPlots, nrow = nrow(tempColDesign), ncol = ncol(tempColDesign), byrow = TRUE)
          if (serpentine) {
               for (j in seq(2, nrow(plotNum[[i]]), by = 2)) { plotNum[[i]][j,] <- rev(plotNum[[i]][j,]) }
          }
          tempNewDesign <- merge(tempNewDesign, as.data.frame.table(plotNum[[i]]), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
          randomize <- data.frame(rbind(randomize, tempNewDesign))
          dimnames(plan[[i]]) <- list(paste("FieldRow",1:numFieldRow, sep = ""),paste("FieldCol",1:numFieldCol, sep = ""))
          dimnames(plotNum[[i]]) <- dimnames(plan[[i]])
     }
     randomize$ID <- as.numeric(as.character(randomize[,"Freq.x"]))
     randomize <- randomize[,c("Trial", "Row", "Column", "ID","Freq.x", "Freq.y", "Var1", "Var2")]

     EntryLabelNew <- NULL
     for (i in (1:trial)) { EntryLabelNew <- c(EntryLabelNew, newTrmt) }
     randomize <- randomize[order(randomize$Trial, randomize$ID),]
     randomize[randomize[,"ID"] %in% 1:numNew, "ID"] <- EntryLabelNew
     for (j in (1:numCheck)) { randomize[randomize[,"ID"] == EntryCodeCheck[j], "ID"] <- checkTrmt[j] }
     
     randomize[,"Var1"] <- as.numeric(randomize[,"Var1"]) # -- fieldRow
     randomize[,"Var2"] <- as.numeric(randomize[,"Var2"]) # -- fieldColumn
     
     if (!is.null(trmtName) && is.valid.name(trimStrings(trmtName[1]))) {
          colnames(randomize) <- c("Trial", "RowBlk", "ColBlk", "ID", trmtName, "PlotNum", "FieldRow", "FieldCol")     
     } else { colnames(randomize) <- c("Trial", "RowBlk", "ColBlk", "ID","EntryNo", "PlotNum", "FieldRow", "FieldCol") }
     
     if (!is.null(trmtLabel) && is.valid.name(trimStrings(trmtLabel[1]))) { 
          colnames(randomize)[match("ID", names(randomize))] <- trmtLabel 
     }
     #randomize <- randomize[order(randomize$Trial,  randomize$FieldRow, randomize$FieldCol), ]
     randomize <- randomize[order(randomize$Trial,  randomize$PlotNum), ]
     
     rownames(randomize) <- 1:nrow(randomize)
     names(plan) <- paste("Trial", 1:trial, sep = "")
     names(plotNum) <- paste("Trial", 1:trial, sep = "")
     
     #if (display) {
          cat(toupper("Design Properties:"),"\n",sep = "")
          cat("\t","Augmented Latin Square Design (Augmented LSD)","\n\n",sep = "")
          cat(toupper("Design Parameters:"),"\n",sep = "")
          cat("\t","Number of Trials = ", trial, "\n",sep = "")
          cat("\t","Number of Replicated Treatments = ", length(checkTrmt), "\n",sep = "")
          cat("\t","Levels of Replicated Treatments = ", sep = "")
          if (length(checkTrmt) <= 5) { cat(paste(checkTrmt, collapse = ", ", sep = ""), "\n", sep = "")
          } else {
               cat(paste(checkTrmt[1:3],collapse = ", ", sep = ""), sep = "")
               cat(", ..., ",checkTrmt[length(checkTrmt)], "\n",sep = "")
          }
          cat("\t","Number of Unreplicated Treatments = ", length(newTrmt), "\n",sep = "")
          cat("\t","Levels of UnReplicated Treatments = ", sep = "")
          if (length(newTrmt) <= 5) { cat(paste(newTrmt, collapse = ", ", sep = ""), "\n", sep = "")
          } else {
               cat(paste(newTrmt[1:3],collapse = ", ", sep = ""), sep = "")
               cat(", ..., ",newTrmt[length(newTrmt)], "\n",sep = "")
          }
          cat("\t","Number of Field Row = ", numFieldRow, "\n",sep = "")
          cat("\t","Number of Field Column = ", numFieldCol, "\n\n",sep = "")
          if (errorDf < 12) {
               cat("WARNING: Too few error df.","\n\n")
          }
          #cat("Results of Randomization:\n")
          #printDataFrame(randomize)
          
     #}
     
     if (!is.null(file)) {
          tempFile <- strsplit(file, split = "\\.")[[1]]
          tempExt <- tolower(tempFile[length(tempFile)])
          if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
          newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
          newFileExt <- paste(newFile, tempExt, sep = ".")
          if (tempExt == "csv") { write.csv(randomize, file = newFileExt, row.names = FALSE)
          } else { save(randomize, file = newFileExt) }
     } else {
          cat("Results of Randomization:\n")
          printDataFrame(randomize)
     }    
     
     # printing of the layout
     for (i in (1:trial)) {
          if (!is.null(file)) { png(filename = paste(newFile, "Layout_Trial",i,".png",sep = ""))  }
          
          if (basicRowWAddlRow == 0 && basicColWAddlCol == 0) {
               des.plot(plan[[i]], seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
                        bwd = 4, bcol = 4, cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
               des.plot(plan[[i]], seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3,
                        bdef = cbind(embedRowInBasicRow, embedColInBasicCol), bwd = 4)
          } else {
               des.plot(plan[[i]], seq(numNew), col = 8, new = TRUE, label = TRUE, chtdiv = 3,
                        bwd = 4, bcol = 4, cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")
               des.plot(plan[[i]], seq(numCheck)+numNew, col = 7, new = FALSE, label = TRUE, chtdiv = 3, bwd = 4)
          }
          if (!is.null(file)) { dev.off() } 
     }
     
     return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum)))
}

     