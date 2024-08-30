# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# printAOVTable
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 08.09.2012
# -------------------------------------------------------------------------------


printAOVTable <- function(aovtable) UseMethod("printAOVTable")

printAOVTable.default <- function(aovtable) {
     
     tempaov <- aovtable
     tempaov[,1] <- format(tempaov[,1], scientific = FALSE, digits = 1)
     tempaov[!is.na(tempaov[,2]), 2] <- format(round(tempaov[!is.na(tempaov[,2]), 2],4), scientific = FALSE, digits = 4, nsmall = 4)
     tempaov[!is.na(tempaov[,3]), 3] <- format(round(tempaov[!is.na(tempaov[,3]), 3],4), scientific = FALSE, digits = 4, nsmall = 4)
     tempaov[(!is.na(tempaov[,4]) | !is.nan(tempaov[,4]) | !is.infinite(tempaov[,4])), 4] <- format(round(tempaov[(!is.na(tempaov[,4]) | !is.nan(tempaov[,4]) | !is.infinite(tempaov[,4])), 4], 2), scientific = FALSE, digits = 2, nsmall = 2)
     tempaov[(!is.na(tempaov[,5]) | !is.nan(tempaov[,5]) | !is.infinite(tempaov[,5])), 5] <- format(round(tempaov[(!is.na(tempaov[,5]) | !is.nan(tempaov[,5]) | !is.infinite(tempaov[,5])), 5], 4), scientific = FALSE, digits = 4, nsmall = 4)
     
	colnames(aovtable) <- c("DF", "Sum of Square", "Mean Square", "F Value", "Pr(> F)")
	colwidth <- max(c(nchar(trimStrings(rownames(aovtable), side = "right")), nchar("Source"))) + 2
     for (i in (1:4)) colwidth <- c(colwidth, max(c(nchar(tempaov[,i]), nchar(colnames(aovtable)[i]))) + 2) 
     colwidth <- c(colwidth, 8)

     options(width = 5000)
	cat(formatC(paste(rep("-", sum(colwidth)), collapse="")), sep = "\n")
	cat(formatC("Source", width = colwidth[1], format = "s", flag = "-"), sep = "")
	for (i in (1:ncol(aovtable))) { cat(formatC(colnames(aovtable)[i], width = colwidth[i + 1], format = "s"), sep = "") }
	cat("\n")
	cat(formatC(paste(rep("-", sum(colwidth)), collapse="")), sep = "\n")
	for (i in (1:nrow(aovtable))) {
	     cat(formatC(trimStrings(rownames(aovtable)[i], side = "right"), width = colwidth[1], format = "s", flag = "-"), sep = "")
	     cat(formatC(tempaov[i,1], width = colwidth[2]), sep = "")
	     for (j in (2:5)) { 
	          if (!is.na(aovtable[i,j])) { cat(formatC(tempaov[i,j], width = colwidth[j + 1], format = "s"), sep = "") 
	          } else { 
                    if (is.nan(aovtable[i,j]) | is.infinite(aovtable[i,j])) {
                         cat(formatC(tempaov[i,j], width = colwidth[j + 1], format = "s"), sep = "")
                    } else {cat(formatC("", width = colwidth[j + 1], format = "s"), sep = "") }
	          }
	     }
	     #for (j in (2:3)) { 
	     #     if (!is.na(aovtable[i,j])) { cat(formatC(aovtable[i,j], digits = 4, width = colwidth[j + 1], format = "f"), sep = "") 
	     #     } else { cat(formatC("", width = colwidth[j + 1], format = "s"), sep = "") }
	     #}
	     #if (!is.na(aovtable[i,4])) { cat(formatC(aovtable[i,4], digits = 2, width = colwidth[5], format = "f"), sep = "") } else { cat(formatC("", width = colwidth[5], format = "s"), sep = "") }
	     #if (!is.na(aovtable[i,5])) { cat(formatC(aovtable[i,5], digits = 4, width = colwidth[6], format = "f"), sep = "") } else { cat(formatC("", width = colwidth[6], format = "s"), sep = "") }
	     
	     cat("\n")
	}
	cat(formatC(paste(rep("-", sum(colwidth)), collapse="")), sep = "\n\n")
} ### end statement