# ----------------------------------------------------------------------
# pairwise.comparison: Function for displaying the mean comparison
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 
# ----------------------------------------------------------------------

pairwiseComparison <- function(aovTable, design, effect, data, respvar, pwTest = NULL, siglevel = 0.05) UseMethod("pairwiseComparison")

pairwiseComparison.default <- function(aovTable, design, effect, data, respvar, pwTest = NULL, siglevel = 0.05) {
  	pwTestAvailable <- c("LSD", "duncan", "HSD", "SNK", "scheffe")
     	#pwLabel <- c("Least Significant Difference (LSD) Test", "Duncan Multiple Range Test (DMRT)", "Tukey's Honest Significant Difference (HSD) Test", "Student-Newmann-Keul's (SNK) Test",  "Scheffe's Test")
     	availableDesign <- c("CRD", "RCBD", "LSD", "SplitCRD", "SplitRCBD", "SplitLSD", "Strip", "Split2CRD", "Split2RCBD", "Split2LSD", "Strip-Split",	"Split3CRD", "Split3RCBD", "Split3LSD", "Strip-Split2")

     	designChoice <- match(design, availableDesign)
     	pwTestChoice <- pwTest
     
     	switch(designChoice,
      	{numfactor <- 1; numblk <- 0}, {numfactor <- 1; numblk <- 1}, {numfactor <- 1; numblk <- 2}, 
            {numfactor <- 2; numblk <- 1}, {numfactor <- 2; numblk <- 1}, {numfactor <- 2; numblk <- 2},
            {numfactor <- 2; numblk <- 1}, {numfactor <- 3; numblk <- 1}, {numfactor <- 3; numblk <- 1},
            {numfactor <- 3; numblk <- 2}, {numfactor <- 3; numblk <- 1}, {numfactor <- 4; numblk <- 1},
            {numfactor <- 4; numblk <- 1}, {numfactor <- 4; numblk <- 2}, {numfactor <- 4; numblk <- 1}
     	)

     	if (length(strsplit(effect, split = ":")[[1]]) == 1) {
      	if (is.null(pwTest)) { if (nlevels(data[,effect]) <= 5) { pwTestChoice <- "LSD" } else { pwTestChoice <- "HSD" }}
          	cat("Pairwise Mean Comparison of ", effect,"\n", sep = "")
          	if (designChoice <= 3) { pairwiseAmong(data, respvar, typeTest = pwTestChoice, trmt = effect, dfError = aovTable[[1]][nrow(aovTable[[1]]),1], MSError = aovTable[[1]][nrow(aovTable[[1]]),3], siglevel) 
          	} else {
            	start <- 1
               	while (is.na(match(effect, trimStrings(rownames(aovTable[[start]][[1]])))) && start <= length(aovTable)) start <- start + 1
		   	pairwiseAmong(data, respvar, typeTest = pwTestChoice, trmt = effect, dfError = aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]),1], MSError = aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]),3], siglevel) 
          	}          
     	} 
     
	# pairwise comparison for interaction effect
     	if (length(strsplit(effect, split = ":")[[1]]) > 1) { 
      	tempFactor <- strsplit(effect, split = ":")[[1]]
          	atLevel <- NULL
          	for (j in (1:length(tempFactor))) { atLevel <- c(atLevel, paste(tempFactor[-I(j)], collapse = ":", sep = "")) }
          	if (designChoice > 3) {
               	tempPos <- NULL
              	for (i in (1:length(tempFactor))) {
                    	cont <- TRUE; index <- 1
                    	while(cont && index <= length(aovTable)){
                         	if (!is.na(match(tempFactor[i], trimStrings(rownames(aovTable[[index]][[1]]))))) { tempPos <- c(tempPos, index) }
                         	index <- index + 1
                    	}
               	}
               	df <- list()
               	MSE <- list()
               	if (length(unique(tempPos)) > 1) {
                    	for (j in (1:length(tempFactor))) {
                         	df[[j]] <- aovTable[[tempPos[j]]][[1]][nrow(aovTable[[tempPos[j]]][[1]]),1]
                         	MSE[[j]] <- aovTable[[tempPos[j]]][[1]][nrow(aovTable[[tempPos[j]]][[1]]),3]
                         	if (any(unique(tempPos[-I(j)]) > tempPos[j])) {
                              	df.num <- aovTable[[tempPos[j]]][[1]][nrow(aovTable[[tempPos[j]]][[1]]),3]
                              	df.denom <- (aovTable[[tempPos[j]]][[1]][nrow(aovTable[[tempPos[j]]][[1]]),3])**2/aovTable[[tempPos[j]]][[1]][nrow(aovTable[[tempPos[j]]][[1]]),1]
                              	prevLevel <- 1
                              	tempEffect <- tempFactor[j]
                              	for (k in (1:length(unique(tempPos[-I(j)])))) {
                                   		if (designChoice == 7 || designChoice == 11 || designChoice == 15) { 
                                        		if (tempPos[j] == 1 || tempPos[j] ==2) { 
                                             		condition <- tempPos[j] != unique(tempPos[-I(j)])[k]
                                             		lowestLevel <- FALSE
                                        		} else {
                                             		condition <- tempPos[j] < unique(tempPos[-I(j)])[k]
                                             		lowestLevel <- TRUE
                                        		}
                                		} else { 
                                        		condition <- tempPos[j] < unique(tempPos[-I(j)])[k]
                                        		lowestLevel <- TRUE 
                                   		}
                                   		if (condition) {
                                        		diffFactor <- tempFactor[-I(j)][which(tempPos[-I(j)] == unique(tempPos[-I(j)])[k])]
                                        		levelTemp <- nlevels(eval(parse(text = paste("data[,'", diffFactor,"']", collapse = ":", sep = ""))))
                                        		tempEffect <- tempFactor[sort(match(c(tempEffect,diffFactor), tempFactor))]
                                        		start <- 1
                                        		while(is.na(match(paste(tempEffect, collapse = ":", sep = ""), trimStrings(row.names(aovTable[[start]][[1]])))) && start <= length(aovTable)) start <- start + 1
                                        		df.num <- df.num + (prevLevel*(levelTemp - 1)*aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]),3])
                                        		df.denom <- df.denom + (((prevLevel*(levelTemp - 1)*aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]),3])**2)/aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]),1])
                                        		prevLevel <- prevLevel * levelTemp
                                   		}   
                             		}
                              	df[[j]] <- df.num**2/df.denom
                              	if (lowestLevel) { MSE[[j]] <- df.num/levelTemp 
                              	} else { MSE[[j]] <- df.num/nlevels(eval(parse(text = paste("data[,'", tempFactor[-I(j)],"']", collapse = ":", sep = "")))) } 
                       		}
               		} ## END FOR STMT -- for (i in (1:length(tempFactor)))
               	} ## END IF STMT -- if (length(unique(tempPos)) > 1)
               	else {
                    	start <- 1
                    	while (is.na(match(effect,trimStrings(row.names(aovTable[[start]][[1]])))) && start <= length(aovTable)) start <- start + 1
                    	df[[1]] <- aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]), 1]
                    	MSE[[1]] <- aovTable[[start]][[1]][nrow(aovTable[[start]][[1]]), 3]
               	} ## END IF ELSE STMT -- if (length(unique(tempPos)) > 1)
      	} ## END IF STMT -- if (designChoice > 3)
          
          	for (j in (1:length(tempFactor))) {
               	if (is.null(pwTest)) { if (nlevels(data[,tempFactor[j]]) <= 5) { pwTestChoice <- "LSD" } else { pwTestChoice <- "HSD" }}
               	cat("Comparison of ",tempFactor[j]," at each level of ", atLevel[j],"\n", sep = "")
               	if (length(strsplit(atLevel[j], split = ":")[[1]]) != 1) {
                    	temp.blevel <- NULL
                    	for (l in (1:length(strsplit(atLevel[j], split = ":")[[1]]))) { temp.blevel[l] <- nlevels(data[,strsplit(atLevel[j], split = ":")[[1]][l]])  }
                    	temp.blevel.max <- strsplit(atLevel[j], split = ":")[[1]][match(max(temp.blevel), temp.blevel)]
                    	temp.blevel.oth <- strsplit(atLevel[j], split = ":")[[1]][-I(match(max(temp.blevel), temp.blevel))]
                    	for (k in (1:length(pwTestChoice))) {
                         	#cat(pwLabel[match(pwTestChoice[k], pwTestAvailable)], "\n\n")
                         	for (l in (1:nlevels(eval(parse(text = paste("data[,'",temp.blevel.oth,"']", collapse = ":", sep = "")))))) {
                              	sub.label1 <- NULL
                              	for (m in (1:length(temp.blevel.oth))) sub.label1 <- c(sub.label1, paste(temp.blevel.oth[m], " = ", strsplit(levels(eval(parse(text = paste("data[,'", temp.blevel.oth,"']", collapse = ":", sep = ""))))[l], split = ":")[[1]][m], sep = ""))
                              	sub.label <- paste(" & (", paste("data[,'",temp.blevel.oth,"']", collapse = ":", sep = ""), ") == levels(", paste("data[,'",temp.blevel.oth,"']", collapse = ":", sep = ""),")[[", l, "]]", sep = "") 
						cat(sub.label1,"\n\n")                              
                              	if (designChoice <= 3) {
							pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,temp.blevel.max]), dfError = aovTable[[1]][nrow(aovTable[[1]]),1], MSError = aovTable[[1]][nrow(aovTable[[1]]),3], f1 = tempFactor[j], f2 = temp.blevel.max, f3 = sub.label, siglevel = siglevel)     
                              	} else {
                              		if (length(unique(tempPos)) > 1) {
								pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,temp.blevel.max]), dfError = df[[j]], MSE[[j]], f1 = tempFactor[j], f2 = temp.blevel.max, f3 = sub.label, siglevel = siglevel)     
                                   		} else {
								pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,temp.blevel.max]), dfError = df[[1]], MSE[[1]], f1 = tempFactor[j], f2 = temp.blevel.max, f3 = sub.label, siglevel = siglevel)     
                                  		}
                              	} 
						cat("\n")
                         	} ### end stmt -- for (l in (1:nlevels(eval(parse(text = paste(temp.blevel.oth, collapse = ":", sep = ""))))))
                    	} ## END FOR STMT -- for (k in (1:length(pwTestChoice)))
               	} else { ## ELSE STMT -- if (length(strsplit()[[1]]) != 1)
             		for (k in (1:length(pwTestChoice))) {
                    		if (designChoice <= 3) {
						pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,atLevel[j]]), dfError = aovTable[[1]][nrow(aovTable[[1]]),1], MSError = aovTable[[1]][nrow(aovTable[[1]]),3], f1 = tempFactor[j], f2 = atLevel[j], siglevel = siglevel)     
                         	} else {
                         		if (length(unique(tempPos)) > 1) {
							pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,atLevel[j]]), dfError = df[[j]], MSError = MSE[[j]], f1 = tempFactor[j], f2 = atLevel[j], siglevel = siglevel)     
                              	} else {
							pairwiseWithin(data, respvar, typeTest = pwTestChoice[k], nobs1 = nlevels(data[,tempFactor[j]]), nobs2 = nlevels(data[,atLevel[j]]), dfError = df[[1]], MSError = MSE[[1]], f1 = tempFactor[j], f2 = atLevel[j], siglevel = siglevel)     
                              	}
                      		}
					cat("\n")
                    	} ## END FOR STMT -- for (k in (1:length(pwTestChoice)))
               	} ## END IF ELSE STMT -- if (length(strsplit()[[1]]) != 1)
          	} ## END FOR STMT -- for (j in (1:length(tempFactor)))
     	} ## END IF STMT -- if (length(strsplit(effect, split = ":")[[1]]) > 1)
}
