# -------------------------------------------------------------
# Statistical Tool for Agricultural Research (STAR)
# -----------------------------------------------------------------------
# printCorrMatrix: Print the correlation Matrix
# Created by: Alaine A. Gulles for International Rice Research Institute
# -----------------------------------------------------------------------

printCorrMatrix <- function(coef, pval = NULL, n) UseMethod("printCorrMatrix")

printCorrMatrix.default <- function(coef, pval = NULL, n) {
	rnameWidth <- max(nchar(rownames(coef))) + 3
	colWidth <- max(nchar(max(n)), max(nchar(rownames(coef))), 6) + 3
	numRow <- 2
	other.label <- NULL
	if (!is.null(pval)) { 
		numRow <- numRow + 1
		other.label <- c(other.label, "p-value") 
	}
	other.label <- c(other.label, "n")
	if (attributes(coef)$class == "table") {
		rname <- NULL
		for (i in (1:ncol(coef))) { 
			rname <- c(rname, paste(formatC(rownames(coef)[i], width = rnameWidth, flag = "-"), formatC("coef", width = 8, flag = "-"), sep = ""), 
				     paste(formatC("", width = rnameWidth, flag = "-"), formatC(rep(other.label), width = 8, flag = "-"), sep = ""))
		}
		cnames <- paste(formatC(rownames(coef), width = colWidth, format = "s"))
		temp.output <- as.table(matrix(NA, nrow = nrow(coef)*numRow, ncol = ncol(coef), dimnames = list(rname, cnames)))
		cnt <- 1
		for (i in (1:ncol(coef))) {
			temp.output[cnt,] <- formatC(coef[i,], digits = 4, width = colWidth, format = "f", drop0trailing = FALSE)
			if (is.na(coef[1,1])) temp.output[cnt,i] <- formatC("", width = colWidth, flag = "-")
			if (!is.null(pval)) {
				for (j in (1:ncol(pval))) {
					if (is.na(pval[i,j])) {	temp.output[cnt+1,j] <- formatC("", format = "s", width = colWidth) }
					else { temp.output[cnt+1,j] <- formatC(pval[i,j], digits = 4, format = "f", width = colWidth, drop0trailing = FALSE) }
				}
				temp.output[cnt+2,] <- formatC(n[i,], width = colWidth)
			} else { temp.output[cnt+1,] <- formatC(n[i,], width = colWidth) }
			cnt <- cnt + numRow
		}
	}
	return(temp.output)
}