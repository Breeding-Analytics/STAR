# -----------------------------------------------------------------------------------
# FactorList: Generate all treatment combination of several factors.
# Created by: Alaine A. Gulles 09.13.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 09.13.2012
# -----------------------------------------------------------------------------------

FactorList <- function(generate) UseMethod("FactorList")

FactorList.default <- function(generate){
	if (!is.list(generate)) stop("The argument 'generate' must be a list.")
	for (i in (1:length(generate))) {
		if (length(generate[[i]]) == 1 && is.numeric(generate[[i]])) {
			generate[[i]] <- paste(names(generate)[i],1:generate[[i]], sep = "")
		} 
	}
	return(generate)
}
