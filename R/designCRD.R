# -------------------------------------------------------------------------------------
# designCRD: Generate randomization for complete block design.
# Created by: Alaine A. Gulles 04.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.08.2013
# -------------------------------------------------------------------------------------

designCRD <- function(generate, r = 2, trial = 1, numFieldRow = 1, serpentine = FALSE, display = TRUE, file = NULL) UseMethod("designCRD")

designCRD.default <- function(generate, r = 2, trial = 1, numFieldRow = 1, serpentine = FALSE, display = TRUE, file = NULL) {

	if (is.null(trial) || trial < 1 || is.character(trial) || length(trial) > 1) { stop("The argument 'trial' should be a single numeric value greater than or equal to 1.") }
	if (is.null(r) || r < 2 || is.character(r) || length(r) > 1) { stop("The argument 'r' should be a single numeric value greater than or equal to 2.") }
	if (missing(generate)) { stop("The argument 'generate' is missing.") }
	if (!is.list(generate)) { stop("The argument 'generate' must be a list.") }
	
	tempComb <- GenerateFactor(generate, times = r)
	if (nrow(tempComb)%%numFieldRow != 0) {   stop("Total number of plot is not divisible by number of field row.")}
	randomize <- NULL
     plan <- list()
     plotNum <- list()
     
     if (numFieldRow == 1) serpentine <- FALSE

     for (i in (1:trial)) {
		temp <- data.frame(Trial = as.character(i), tempComb, PlotNum = sample(nrow(tempComb), nrow(tempComb), replace = FALSE))
		temp <- temp[order(temp[,"PlotNum"]),]
		if (ncol(tempComb) > 1) { trmtLabel <- eval(parse(text = paste("paste(temp[,'", paste(names(temp)[2:(ncol(temp)-1)], collapse = "'],' ',temp[,'", sep = ""),"'], sep = '')", sep = "")))
		} else { trmtLabel <- temp[,2] }
          plan[[i]] <- matrix(trmtLabel, nrow = numFieldRow, ncol = nrow(tempComb)/numFieldRow, byrow = TRUE)
		plotNum[[i]] <- matrix(temp[,ncol(temp)], nrow = numFieldRow, ncol = nrow(tempComb)/numFieldRow, byrow = TRUE)
          if (serpentine) {
               for (j in seq(2, numFieldRow, by = 2)) { 
                    plotNum[[i]][j,] <- rev(plotNum[[i]][j,]) 
                    plan[[i]][j,] <- rev(plan[[i]][j,]) 
               }
          }
          tempOrder <- as.data.frame.table(plotNum[[i]])
          tempOrder <- tempOrder[order(tempOrder[,"Freq"]),]
          temp <- cbind(temp, as.numeric(tempOrder[,1]), as.numeric(tempOrder[,2]))
          colnames(temp)[(ncol(temp)-1):ncol(temp)] <- c("FieldRow", "FieldColumn")
		randomize <- rbind(randomize, temp)
	}
     randomize <- randomize[order(randomize$Trial, randomize$PlotNum),]
	rownames(randomize) <- 1:nrow(randomize)
     names(plan) <- paste("Trial", 1:trial, sep = "")
     names(plotNum) <- paste("Trial", 1:trial, sep = "")

	## display in the console the output
	if (display) { 
		cat(toupper("Design Properties:"),"\n",sep = "")
		if (ncol(tempComb) == 1) { cat("\t","Single Factor","\n",sep = "") } else { cat("\t","Factorial Design","\n",sep = "") }
		cat("\t","Completely Randomized Design","\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Replicates = ", r, "\n",sep = "")
		if (ncol(tempComb) == 1) {
			cat("\t","Treatment Name = ", names(tempComb)[1], "\n",sep = "")
			cat("\t","Treatment Levels = ", sep = "")
			if (nlevels(tempComb[,1]) <= 5) { cat(paste(levels(tempComb[,1]), collapse = ", ", sep = ""), sep = "")
			} else {
				cat(paste(levels(tempComb[,1])[1:3], collapse = ", ", sep = ""), sep = "")
				cat(paste(", ...,", levels(tempComb[,1])[nlevels(tempComb[,1])]), sep = "")
			}
			cat("\n\n")
		} else {
			for (i in (1:ncol(tempComb))) {
				cat("\t","Factor ",i," = ", names(tempComb)[i], "\n",sep = "")
				cat("\t","Levels = ", sep = "")
				if (nlevels(tempComb[,i]) <= 5) { cat(paste(levels(tempComb[,i]), collapse = ", ", sep = ""), sep = "")
				} else {
					cat(paste(levels(tempComb[,i])[1:3], collapse = ", ", sep = ""), sep = "")
					cat(paste(", ...,", levels(tempComb[,i])[nlevels(tempComb[,i])]), sep = "")
				}
				cat("\n")
			}
			cat("\n")
		}
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
