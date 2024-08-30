# -------------------------------------------------------
# AppendData: Function for appending two data frames
# Created by: Alaine A. Gulles 04.11.2012 for
#	      International Rice Research Institute
# Modified by: Alaine A. Gulles 08.28.2012
# -------------------------------------------------------

AppendData <- function(MasterData, TransactionData, byMaster, byTransact, MasterVarKeep = NULL, TransactVarKeep = NULL) UseMethod("AppendData")

AppendData.default <- function(MasterData, TransactionData, byMaster, byTransact, MasterVarKeep = NULL, TransactVarKeep = NULL) {

	if (is.character(MasterData)) { MasterData <- eval(parse(text = MasterData)) }
	if (is.character(TransactionData)) { TransactionData <- eval(parse(text = TransactionData)) }

	MasterData <- MasterData[,c(byMaster, MasterVarKeep)]
	TransactionData <- TransactionData[,c(byTransact, TransactVarKeep)]
	names(TransactionData)[1:length(byTransact)] <- names(MasterData)[1:length(byMaster)] 

	if (!is.null(MasterVarKeep) && !is.null(TransactVarKeep)) {
		if (length(intersect(MasterVarKeep,TransactVarKeep)) != 0) {
			commonVar <- intersect(MasterVarKeep,TransactVarKeep)
			index <- which(commonVar == MasterVarKeep)
			MasterVarKeep[index] <- paste(MasterVarKeep[index], ".1", sep = "")
			names(MasterData)[which(commonVar == names(MasterData))] <- MasterVarKeep[index]
			index <- which(commonVar == TransactVarKeep)
			TransactVarKeep[index] <- paste(TransactVarKeep[index], ".2", sep = "")
			names(TransactionData)[which(commonVar == names(TransactionData))] <- TransactVarKeep[index]
			remove(list = c("index", "commonVar"))
		}
	}

	if (!is.null(MasterVarKeep)) {
		if (!is.null(TransactVarKeep)) {
			TransactionData <- cbind(TransactionData[c(1:length(byTransact))], 
				 		 	 matrix(data = NA, nrow = nrow(TransactionData), ncol = length(MasterVarKeep), dimnames = list(NULL, MasterVarKeep)),
						 	 TransactionData[,c((length(byTransact)+1):ncol(TransactionData))])
			
			names(TransactionData) <- c(byMaster, MasterVarKeep, TransactVarKeep)
		} else {
			TransactionData <- cbind(TransactionData[,c(1:length(byTransact))], 
				 		 	 matrix(data = NA, nrow = nrow(TransactionData), ncol = length(MasterVarKeep), dimnames = list(NULL, MasterVarKeep)))

 		}
	} 

	if (!is.null(TransactVarKeep)) {
		MasterData <- cbind(MasterData, matrix(data = NA, nrow = nrow(MasterData), ncol = length(TransactVarKeep), dimnames = list(NULL, TransactVarKeep)))
	}
	appendData <- rbind(data.frame(MasterData, Source = 1), data.frame(TransactionData, Source = 2))
	rownames(appendData) <- 1:nrow(appendData)
	return(appendData)
}
