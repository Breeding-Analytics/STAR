# --------------------------------------------------------------
# Created by AAGulles for International Rice Research Institute
# --------------------------------------------------------------

printTable <- function(table) UseMethod("printTable")

printTable.default <- function(table) {
	if(!is.table(table)) { stop("The argument 'table' should be a class of table.") }
	#if (length(dimnames(table)) != 2)) { stop("The function works only for n X r") }

	colwidth1 <- max(nchar(names(dimnames(table)[1])), nchar(dimnames(table)[[1]])) + 4
	colwidth2 <- max(max(nchar(table), nchar(dimnames(table)[[2]])), nchar(dimnames(table)[[2]])) + 2

	tableBorder <- function(withLine = TRUE) {
		if (withLine) { cat(paste(rep("-", colwidth1), collapse = ""),"+",sep = "")
		} else { cat(formatC(" ", width = colwidth1, format = "s", flag = "-"), "|",sep = "") }
		for (i in (1:length(dimnames(table)[[2]]))) {
			cat(paste(rep("-", colwidth2), collapse = ""),"-",sep = "")
		}
		cat("\n")
	}

	tableBorder(withLine = TRUE)
	
	cat(formatC(" ", width = colwidth1, format = "s", flag = "-"), "|",sep = "")
	for (i in (1:length(dimnames(table)[[2]]))) {
		if (i == 1) { cat(formatC(names(dimnames(table)[2]), width = colwidth2+1, format = "s", flag = "-"),sep = "")
		} else { 
			if (i == length(dimnames(table)[[2]])) { cat(formatC(" ", width = colwidth2, format = "s", flag = "-"),"|",sep = "")
			} else { cat(formatC(" ", width = colwidth2+1, format = "s", flag = "-"),sep = "")	}
		} 
	}

	cat("\n")

	tableBorder(withLine = FALSE)

	cat(formatC(names(dimnames(table)[1]), width = colwidth1, format = "s", flag = "-"), "|",sep = "")
	for (i in (1:length(dimnames(table)[[2]]))) {
		cat(formatC(colnames(table)[i], width = colwidth2, format = "s", flag = "-"),"|",sep = "")
	}
	cat("\n")

	tableBorder(withLine = TRUE)
	for (j in (1:nrow(table))) {
		cat(formatC(rownames(table)[j],width = colwidth1, format = "s", flag = "-"),"|",sep = "")
		for (i in (1:length(dimnames(table)[[2]]))) {
			cat(formatC(table[j,i], width = colwidth2, format = "d"),"|",sep = "")
		}
		cat("\n")
	}

	tableBorder(withLine = TRUE)
}
