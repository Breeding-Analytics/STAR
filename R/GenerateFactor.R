# -----------------------------------------------------------------------------------
# GenerateFactor: Generate all treatment combination of several factors.
# Created by: Alaine A. Gulles 05.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.11.2012
# -----------------------------------------------------------------------------------

GenerateFactor <- function(generate, times = 1) UseMethod("GenerateFactor")

GenerateFactor.default <- function(generate, times = 1) {
	if (!is.list(generate)) stop("The argument 'generate' must be a list.")	
	numFactor <- length(generate)
	trmtNumLevel <- NULL

	generate <- FactorList(generate)

	if (length(generate) == 1) { trmt <- data.frame(gl(length(generate[[1]]), times, labels = generate[[1]])) 
	} else {
		for (i in (1:length(generate))) { trmtNumLevel <- c(trmtNumLevel, length(generate[[i]])) }
		trmt <- data.frame(gl(length(generate[[1]]), prod(trmtNumLevel[2:length(trmtNumLevel)])*times, labels = generate[[1]]))
		for (i in (2:length(trmtNumLevel))) {
			if (i == length(trmtNumLevel)) { trmt <- cbind(trmt, data.frame(gl(trmtNumLevel[i], times, prod(trmtNumLevel[1:length(trmtNumLevel)])*times, labels = generate[[i]])))
			} else { trmt <- cbind(trmt, data.frame(gl(trmtNumLevel[i], prod(trmtNumLevel[(i+1):length(trmtNumLevel)])*times, prod(trmtNumLevel[1:length(trmtNumLevel)])*times, labels = generate[[i]]))) }
		} 
	} 
	
	names(trmt) <- names(generate)
	return(trmt)
}
