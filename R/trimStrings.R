# -------------------------------------------------------------------------------
# RCropStat Utilities: Function
# -------------------------------------------------------------------------------
# trimStrings: Trim whitespace from start and/or end of the string
# Adapted from the function trim.blanks of Rcmdr package
# Modified by: Alaine A. Gulles for International Rice Research Institute
# -------------------------------------------------------------------------------
# Agruments: text - character vector
#            side - side on which whitespace will be removed
#                 - possible values: c("both", "left", "right")
# -------------------------------------------------------------------------------

trimStrings <- function(text, side = "both") {
	if (length(side) != 1) side <- side[1]
	if (is.na(match(side, c("both", "left", "right")))) { side <- "both" }
	if (side == "right") { gsub("\ *$", "", text)	}
	else { if (side == "left") { gsub("^\ *", "", text) }
		 else { gsub("^\ *", "", gsub("\ *$", "", text)) }
	}
}

