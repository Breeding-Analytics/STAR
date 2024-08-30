# --------------------------------------------------------
# RCropStat Beta Menu
# Functions obtained from package Rcmdr
# --------------------------------------------------------

is.valid.name <- function(x) { 
	length(x) == 1 && is.character(x) && x == make.names(x) 
}
