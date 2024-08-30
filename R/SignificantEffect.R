# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# SignificantEffect: determine the significant effect from the ANOVA table
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.19.2013
# -------------------------------------------------------------------------------

SignificantEffect <- function(aovTable, alpha = 0.05) UseMethod("SignificantEffect")

SignificantEffect.default <- function(aovTable, alpha = 0.05) {
     rownames(aovTable) <- trimStrings(rownames(aovTable))  # remove spaces of sv in anova table
     tableSigEffect <- aovTable[!is.na(aovTable[,5]),]
     tempEffect <- NULL
     if (all(tableSigEffect[,5] > alpha)) { sigEffect <- NULL
     } else {
          sigEffect <- NULL
          hEffect <- rownames(tableSigEffect)[nrow(tableSigEffect)]
          allFactor <- strsplit(hEffect, split = ":")[[1]]
          tableSigEffect <- tableSigEffect[tableSigEffect[,5] <= alpha,]
          numOfEffect <- 0; Freq <- 0
          if (!is.na(match(hEffect, rownames(tableSigEffect)))) { sigEffect <- hEffect 
          } else {
               if (length(allFactor) > 1) { 
                    tableSigEffect <- cbind(tableSigEffect, numOfEffect = unlist(lapply(strsplit(rownames(tableSigEffect),split = ":"),"length")))
                    numFactorCheck <- max(tableSigEffect[,"numOfEffect"])
                    if (numFactorCheck > 1) { checkOtherEffect <- TRUE 
                    } else { 
                         checkOtherEffect <- FALSE 
                         #tableSigEffect <- tableSigEffect[-I(match(trimStrings(rownames(tableSigEffect))[-I(match(trimStrings(allFactor), trimStrings(rownames(tableSigEffect))))], trimStrings(rownames(tableSigEffect)))),]
                    }
                    # determine intial significant effect
                    sigEffect <- c(sigEffect, trimStrings(rownames(subset(tableSigEffect, numOfEffect == numFactorCheck))))
                    tempEffect <- strsplit(paste(sigEffect, collapse = ":", sep = ""), split = ":")[[1]]
                    tempEffect <- data.frame(Freq = tapply(tempEffect, tempEffect, length))
                    while (checkOtherEffect) {
                         factorToCheck <- NULL
                         if (length(unique(rownames(tempEffect))) < length(allFactor)) {
                              tempFactor <- allFactor[-I(sort(match(rownames(tempEffect),allFactor)))]
                         } else {
                              if (nrow(subset(tempEffect, Freq == 1)) > 0) {
                                   tempFactor <- allFactor[sort(match(rownames(subset(tempEffect, Freq == 1)), allFactor))]
                              } else { checkOtherEffect <- FALSE; tempFactor <- NULL }
                         }
                         if (length(tempFactor) > 1) {
                              numFactorCheck <- numFactorCheck - 1
                              if (numFactorCheck <= 1) { checkOtherEffect <- FALSE 
                              } else {
                                   if (length(tempFactor) < numFactorCheck) numFactorCheck <- length(tempFactor)
                                   for (i in (1:ncol(combn(tempFactor,numFactorCheck)))) {
                                        factorToCheck <- c(factorToCheck,paste(combn(tempFactor, numFactorCheck)[,i], collapse = ":", sep = ""))
                                   }
                              }
                         } else { factorToCheck <- tempFactor; checkOtherEffect <- FALSE }
                         
                         if (!is.null(factorToCheck)) {
                              newSigFactor <- FALSE
                              for (i in (1:length(factorToCheck))) {
                                   if (nrow(subset(tableSigEffect, trimStrings(rownames(tableSigEffect)) == factorToCheck[i])) != 0) {
                                        sigEffect <- c(sigEffect,trimStrings(rownames(subset(tableSigEffect, trimStrings(rownames(tableSigEffect)) == factorToCheck[i]))))
                                        newSigFactor <- TRUE
                                   }
                              }
                              if (newSigFactor) {
                                   tempEffect <- strsplit(paste(sigEffect, collapse = ":", sep = ""), split = ":")[[1]]
                                   tempEffect <- data.frame(Freq = tapply(tempEffect, tempEffect, length))
                              }
                         } ## END IF STMT -- if (!is.null(factorToCheck))
                    } ## END WHILE STMT
               } ## END IF STMT -- if (length(allFactor) > 1)
          } ## END IF ELSE STMT -- if (!is.na(match(hEffect, rownames(tableSigEffect))))
     } ## END STMT -- if (all(tableSigEffect[,5] > alpha))
     return(sigEffect)
} ## END FUNCTION
