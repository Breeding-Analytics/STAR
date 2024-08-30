# -------------------------------------------------------------------------------------
# designPRep
# Description: Generate randomization and layout for augmented Row-Column design.
# Script Created by: Violeta I Bartolome for International Rice Research Institute
# Function Created by: Alaine A. Gulles 10.15.2014 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.22.2015
# Note: Uses the DiGGer function from package DiGGer
# -------------------------------------------------------------------------------------
# Arguments:
# numCheck
# -------------------------------------------------------------------------------------

designPRep <- function(trmtPerGrp, trmtRepPerGrp, trmtName = NULL, trial = 1, numFieldRow, 
                       serpentine = FALSE, trmtLabel = NULL, trmtListPerGrp = NULL, file = NULL) UseMethod("designPRep")

designPRep.default <- function(trmtPerGrp, trmtRepPerGrp, trmtName = NULL, trial = 1, numFieldRow, 
                               serpentine = FALSE, trmtLabel = NULL, trmtListPerGrp = NULL, file = NULL) {
     
     # --- check input from user:

     if (length(trmtPerGrp) < 2) { stop("Error: At least two groups should be entered.") }
     if (length(trmtPerGrp) != length(trmtRepPerGrp)) { stop("Error: Number of groups should be of the same length as the length of number of replicates per group.") }
     numTrmtPerGrp <- unlist(trmtPerGrp)
     totNumTrmt <- sum(numTrmtPerGrp)
     
     #if (blk < 2) {  stop("The number of blocks should be at most 2.") }
     
     if (totNumTrmt > 1500) {  stop("The maximum number of experimental units that can be generated is 1500.") }
     
     #if (max(trmtRepPerGrp) > blk) { stop("Error: Max number of replicates per group should be at most the number of blocks.") }
     
     trmt <- NULL
     #numTrmtPerGrp <- NULL
     if (is.null(trmtListPerGrp)){
          for (i in (1:length(trmtPerGrp))) {
               trmt <- rbind(trmt, data.frame(GRP = names(trmtPerGrp)[i],
                                          TRMTLEVEL = paste(names(trmtPerGrp)[i], 1:trmtPerGrp[[i]], sep = "_")))
               #numTrmtPerGrp <- c(numTrmtPerGrp, trmtPerGrp[[i]])
          }
          
     } else {
          if (length(trmtListPerGrp) != totNumTrmt) { stop("Error: ") }
          trmt <- data.frame(GRP = rep(names(trmtPerGrp), times = numTrmtPerGrp), TRMTLEVEL = trmtListPerGrp)
     }
     
     trmt[,"TRMTLEVEL"] <- as.character(trmt[,"TRMTLEVEL"])
     totPlots <- as.numeric(crossprod(numTrmtPerGrp,trmtRepPerGrp))
     numFieldCol <- totPlots/numFieldRow
     if (!is.wholenumber(numFieldCol)) { stop("Error: Total number of plots should be divisible by the number of field rows.") }

     #numPlotPerBlk <- totPlots/blk
     #if (!is.wholenumber(numPlotPerBlk)) { stop("Error: Total number of plots should be divisible by the number of blocks.") }
     #colPerBlk <- numPlotPerBlk/rowPerBlk
     #if (!is.wholenumber(colPerBlk)) { stop("Error: Total number of plots should be divisible by the number of blocks.") }
     #numBlkRow <- numFieldRow/rowPerBlk
     #if (!is.wholenumber(numBlkRow)) { stop("Error: Total number of plots should be divisible by the number of row per blocks.") }
     #numBlkCol <- numFieldCol/colPerBlk
     #if (!is.wholenumber(numBlkCol)) { stop("Error: Total number of plots should be divisible by the number of row per blocks.") }
     
     randomize <- NULL
     plan <- list()
     plan1 <- list()
     plotNum <- NULL
     #blkNum <- NULL
          
     if(length(numTrmtPerGrp) <= 8) { colorCode <- c(8:(8-length(numTrmtPerGrp) + 1))     
     } else { colorCode <- rev(c(1:length(numTrmtPerGrp))) } 
     
     if (!is.null(file)) { 
          tempFile <- strsplit(file, split = "\\.")[[1]]
          tempExt <- tolower(tempFile[length(tempFile)])
          if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
          newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
          #newFile <- paste(newFile, tempExt, sep = ".")          
     }
     
     for (i in (1:trial)) {
          
          result <- try(prep <- DiGGer(NumberOfTreatments = totNumTrmt,
                                       RowsInDesign = numFieldRow,
                                       ColumnsInDesign = numFieldCol,
                                       #BlockIn2D = c(rowPerBlk, colPerBlk),           
                                       TreatmentRepeatsPerReplicate = rep(trmtRepPerGrp, numTrmtPerGrp),
                                       TreatmentName = trmt[,"TRMTLEVEL"],
                                       TreatmentNumber = c(1:totNumTrmt)), 
                        silent = TRUE)
          
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }
          
          capture.output(run(prep))
          if(all(class(result) == "try-error")) {
               msg <- trimStrings(strsplit(result, ":")[[1]])
               msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
               stop(paste("Error in DiGGer:", msg, sep = ""))
          }

          prep_mat <- getDesign(prep)

          # --- creating the field layout
          if (!is.null(file)) { png(filename = paste(newFile, "Layout_Trial",i,".png",sep = ""))  }
          for (j in (1:length(trmtPerGrp))) {
               if (j == 1) { des.plot(prep_mat, seq(numTrmtPerGrp[j]), col = colorCode[j], new = TRUE, chtdiv = 3,
                                      cstr = paste("Layout for Trial ",i,":\n\nFieldCol", sep = ""), rstr = "FieldRow")     
               } else {
                    if ( j == length(trmtPerGrp)) { 
                         des.plot(prep_mat, seq(numTrmtPerGrp[j]) + sum(numTrmtPerGrp[1:(j-1)]), col = colorCode[j], 
                                  #new = FALSE, chtdiv = 3, bdef = cbind(rowPerBlk, colPerBlk), bwd = 4)     
                                  new = FALSE, chtdiv = 3, bwd = 4)     
                    } else { des.plot(prep_mat, seq(numTrmtPerGrp[j]) + sum(numTrmtPerGrp[1:(j-1)]), 
                                      col = colorCode[j], new = FALSE, chtdiv = 3)          
                    }
               }
          } ## end stmt -- for (j in (1:length(trmtPerGrp)))
          if (!is.null(file)) { dev.off() } 
          
          
          plan[[i]] <- matrix(print(prep, option = "list")$ENTRY, nrow(prep_mat), ncol(prep_mat))
          plan1[[i]] <- matrix(print(prep, option = "list")$ID, nrow(prep_mat), ncol(prep_mat))
          
          if (i == 1) {
               #tempPlotNumBlk <- matrix(1:numPlotPerBlk, rowPerBlk, colPerBlk, byrow = TRUE)
               #if (serpentine) { for (k in seq(2, rowPerBlk, by = 2)) { tempPlotNumBlk[k,] <- rev(tempPlotNumBlk[k,]) }}
               #blkNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               #plotNum <- matrix(0, nrow(plan[[i]]), ncol(plan[[i]]))
               
               #blkcode <- 1
               #for (j in 1:numBlkRow) {
               #     rowIndexLL <- (j * rowPerBlk) - rowPerBlk + 1
               #     rowIndexUL <- rowIndexLL + rowPerBlk - 1
               #     for (k in 1:numBlkCol) {
               #          colIndexLL <- (k * colPerBlk) - colPerBlk + 1
               #          colIndexUL <- colIndexLL + colPerBlk - 1
               #          
               #          blkNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- blkcode
               #          plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- as.numeric(paste(blkcode,paste(c(rep(0, nchar(numPlotPerBlk))), collapse = ""), sep = "")) + tempPlotNumBlk
               #          blkcode <- blkcode + 1 
               #     } ## end stmt -- for (k in 1:numRepCol)
               #} ## end stmt -- for (j in 1:numRepRow)     
               
               plotNum <- matrix(1:totPlots, nrow(plan[[i]]), ncol(plan[[i]]), byrow = TRUE)
               if (serpentine) { for(k in seq(2, numFieldRow, by = 2)) { plotNum[k,] <- rev(plotNum[k,]) }}
               
          } ## end stmt -- if (i == 1)
          
          #tempFieldOrder <- merge(as.data.frame.table(blkNum), as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), suffixes = c("BLK", "PLOTNO")) 
          tempFieldOrder <- as.data.frame.table(plotNum)
          tempFieldOrder[,"Var1"] <- as.numeric(tempFieldOrder[,"Var1"])
          tempFieldOrder[,"Var2"] <- as.numeric(tempFieldOrder[,"Var2"])
          randomize <- rbind(randomize, cbind(TRIAL = i, merge(print(prep, option = "list"), tempFieldOrder, by.x = c("ROW", "RANGE"), by.y = c("Var1", "Var2"))))
          dimnames(plan[[i]]) <- list(paste("FIELDROW", 1:nrow(plan[[i]]), sep = ""), paste("FIELDCOL", 1:ncol(plan[[i]]), sep = ""))
          dimnames(plan1[[i]]) <- dimnames(plan[[i]])

     } ## end stmt -- for (i in (1:trial))
     
     dimnames(plotNum) <- dimnames(plan[[1]])
     names(plan) <- paste("Trial", 1:trial, sep = "")
     names(plan1) <- names(plan)
     randomize <- randomize[,c("TRIAL", "ID", "ENTRY", "Freq", "ROW", "RANGE")]
     #colnames(randomize) <- gsub("Freq","", colnames(randomize))
     colnames(randomize)[match("TRIAL", names(randomize))] <- "Trial"
     colnames(randomize)[match("Freq", names(randomize))] <- "PlotNum"
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
     cat("\t", "p-Rep Design", "\n\n", sep = "")
     cat(toupper("Design Parameters"), "\n", sep = "")
     cat("\t", "Number of Trials = ", trial, "\n", sep = "")
     cat("\t", "Number of Groupings = ", length(trmtPerGrp), "\n", sep = "")
     for (i in (1:length(trmtPerGrp))) {
          cat("\t", "Group ", i, ": ", names(trmtPerGrp)[i], "\n", sep = "")
          cat("\t", "Number of Replicate for ", names(trmtPerGrp)[i], " = ",  trmtRepPerGrp[i],"\n", sep = "")
          cat("\t", "Levels of ", names(trmtPerGrp)[i]," = ", sep = "")
          if (!is.null(trmtListPerGrp)) {
               if (i == 1) {
                    start <- 1 
                    end <- numTrmtPerGrp[[1]]
               } else {
                    start <- sum(numTrmtPerGrp[1:(i-1)])+1
                    end <- start + numTrmtPerGrp[i]-1
               }
               if (trmtPerGrp[[i]] <= 5) { 
                    cat(paste(trmtListPerGrp[start:end], collapse = ", ", sep = ""), "\n", sep = "") 
               } else {
                    cat(paste(trmtListPerGrp[start:(start+2)], collapse = ", ", sep = ""), sep = "") 
                    cat(", ..., ", paste(trmtListPerGrp[end], sep = ""), "\n", sep = "") 
               }
          } else {
               if (trmtPerGrp[[i]] <= 5) { cat(paste(names(trmtPerGrp)[i], paste(1:trmtPerGrp[[i]]), collapse = ", ", sep = ""), "\n", sep = "") 
               } else {
                    cat(paste(names(trmtPerGrp)[i], paste(1:3), collapse = ", ", sep = ""), sep = "") 
                    cat(", ..., ", paste(names(trmtPerGrp)[i], trmtPerGrp[[i]], sep = ""), "\n", sep = "") 
               }
          }
          cat("\n")
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
     
     if (is.null(trmtLabel)) { return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum)))
     } else { return(invisible(list(fieldbook = randomize, layout = plan, layoutID = plan1, plotNum = plotNum))) }
     
} ## end stmt -- designPRep