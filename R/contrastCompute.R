# -------------------------------------------------------------------------------
# R-CropStat Beta Version: Functions for ANALYZE - ANALYSIS OF VARIANCE SUBMENU
# -------------------------------------------------------------------------------
# ContrastCompute: Function for partitioning sum of squares
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 
# -------------------------------------------------------------------------------

#data = tempData
#aovTable = tempAnova
#contrastOption <- list(factor = c("das","depth"), type = c("control","control"), level = c("14","5"))
#contrastOption <- list(factor = c("das","depth"), type = c("control","control"), level = c("14","5"))

ContrastCompute <- function(data, aovTable, mymodel, mymodel2, contrast.option) UseMethod("ConstrastCompute")

ContrastCompute.default <- function(data, aovTable, mymodel, mymodel2, contrast.option) {
	
	cm <- list()
	cm.label1 <- NULL
	cm.label2 <- NULL
	for (i in (1:length(contrast.option$factor))) {
		if (contrast.option$type[i]  == "control") {
			cm[[i]] <- t(contr.treatment(nlevels(eval(parse(text = paste("data$", contrast.option$factor[i], sep = "")))), base = match(contrast.option$level[i], levels(eval(parse(text = paste("data$", contrast.option$factor[i], sep = "")))))))
			for (j in (1:nrow(cm[[i]]))) { rownames(cm[[i]])[j] <- paste(contrast.option$level[i], "vs", levels(eval(parse(text = paste("data$", contrast.option$factor[i], sep = ""))))[-c(match(contrast.option$level[i], levels(eval(parse(text = paste("data$", contrast.option$factor[i], sep = ""))))))][j]) }
		} else {
			if (contrast.option$type[i] == "orthoPoly") {
				cm[[i]] <- t(poly(as.numeric(levels(eval(parse(text = paste("data$", contrast.option$factor[i], sep = ""))))), degree = as.numeric(contrast.option$level[i])))
				#for (j in (1:nrow(cm[[i]]))) { rownames(cm[[i]])[j] <- paste(contrast.option$level[i], "vs", levels(eval(parse(text = contrast.option$factor[i])))[-c(match(contrast.option$level[i], levels(eval(parse(text = contrast.option$factor[i])))))][j]) }
			} else {
				if (contrast.option$type[i] == "user") {
					if (!is.null(contrast.option$coef[[i]])) cm[[i]] <- contrast.option$coef[[i]]
				}
			}
		}
		
		cm.label1[i] <- paste(contrast.option$factor[i], " = make.contrasts(cm[[",i,"]])", sep = "")
		cm.label.temp2 <- NULL
		for (j in (1:nrow(cm[[i]]))) { cm.label.temp2 <- paste(c(cm.label.temp2, paste("'", rownames(cm[[i]])[j], "' = ", j, sep = "")), collapse = " ,", sep = "") }
		cm.label2[i] <- paste(contrast.option$factor[i], " = list(", cm.label.temp2, ")", sep = "")
		
	}

	result.aovcontrast <- suppressWarnings(eval(parse(text = paste("aov(", mymodel, ", contrasts = list(", paste(cm.label1, collapse = ",", sep = ""), "), data = data)", sep = ""))))
	summary.aovcontrast <- eval(parse(text = paste("summary(result.aovcontrast, split = list(", paste(cm.label2, collapse = ",", sep = ""), "))", sep = "")))
     
     if (attr(summary.aovcontrast,"class")[[1]] == "summary.aovlist") {
          tempResult <- eval(parse(text = paste("aov(", mymodel2, ", contrasts = list(", paste(cm.label1, collapse = ",", sep = ""), "), data = data)", sep = "")))
          tempSummary <- eval(parse(text = paste("summary(tempResult, split = list(", paste(cm.label2, collapse = ",", sep = ""), "))", sep = "")))
          for(i in (2:length(summary.aovcontrast))) {
               for (j in (1:2)) {
                    if (j == 1) { condition <- trimStrings(rownames(summary.aovcontrast[[i]][[1]][is.na(summary.aovcontrast[[i]][[1]][,"Sum Sq"]),]), side = "right")
                    } else { condition <- trimStrings(rownames(summary.aovcontrast[[i]][[1]][summary.aovcontrast[[i]][[1]][,"Sum Sq"] == 0,]), side = "right") }
                    if (length(condition) > 0) {
                         index1 <- match(condition, trimStrings(rownames(summary.aovcontrast[[i]][[1]]), side = "right"))
                         index2 <- match(condition, trimStrings(rownames(tempSummary[[1]]), side = "right"))  
                         summary.aovcontrast[[i]][[1]][index1,] <- tempSummary[[1]][index2,]     
                    }
               }
               if (summary.aovcontrast[[i]][[1]][nrow(summary.aovcontrast[[i]][[1]]),1] != aovTable[[i]][[1]][nrow(aovTable[[i]][[1]]),1]) {
                    summary.aovcontrast[[i]][[1]][nrow(summary.aovcontrast[[i]][[1]]),1] <- aovTable[[i]][[1]][nrow(aovTable[[i]][[1]]),1]
                    summary.aovcontrast[[i]][[1]]["Mean Sq"] <- summary.aovcontrast[[i]][[1]]["Sum Sq"]/summary.aovcontrast[[i]][[1]]["Df"]
                    summary.aovcontrast[[i]][[1]][1:(nrow(summary.aovcontrast[[i]][[1]])-1),"F value"] <- summary.aovcontrast[[i]][[1]][1:(nrow(summary.aovcontrast[[i]][[1]])-1),"Mean Sq"]/summary.aovcontrast[[i]][[1]][nrow(summary.aovcontrast[[i]][[1]]),"Mean Sq"]
                    summary.aovcontrast[[i]][[1]][1:(nrow(summary.aovcontrast[[i]][[1]])-1),"Pr(>F)"] <- pf(summary.aovcontrast[[i]][[1]][1:(nrow(summary.aovcontrast[[i]][[1]])-1),"F value"], summary.aovcontrast[[i]][[1]][1:(nrow(summary.aovcontrast[[i]][[1]])-1),"Df"], aovTable[[i]][[1]][nrow(aovTable[[i]][[1]]),1], lower.tail = FALSE)
               }
          }
     } else {
          if (summary.aovcontrast[[1]][nrow(summary.aovcontrast[[1]]),"Df"] != aovTable[[1]][nrow(aovTable[[1]]),"Df"]) {
               summary.aovcontrast[[1]][nrow(summary.aovcontrast[[1]]),"Df"] <- aovTable[[1]][nrow(aovTable[[1]]),"Df"]
               summary.aovcontrast[[1]]["Mean Sq"] <- summary.aovcontrast[[1]]["Sum Sq"]/summary.aovcontrast[[1]]["Df"]
               summary.aovcontrast[[1]][1:(nrow(summary.aovcontrast[[1]])-1),"F value"] <- summary.aovcontrast[[1]][1:(nrow(summary.aovcontrast[[1]])-1),"Mean Sq"]/summary.aovcontrast[[1]][nrow(summary.aovcontrast[[1]]),"Mean Sq"]
               summary.aovcontrast[[1]][1:(nrow(summary.aovcontrast[[1]])-1),"Pr(>F)"] <- pf(summary.aovcontrast[[1]][1:(nrow(summary.aovcontrast[[1]])-1),"F value"], summary.aovcontrast[[1]][1:(nrow(summary.aovcontrast[[1]])-1),"Df"], aovTable[[1]][nrow(aovTable[[1]]),1], lower.tail = FALSE)
          }
     } 
	AOVContrastTable <- ConstructAOVTable(summary.aovcontrast)
	AOVContrastTable[nrow(AOVContrastTable),] <- ConstructAOVTable(aovTable)[nrow(ConstructAOVTable(aovTable)),]
     
	# --- PRINT ANOVA TABLE --- #
	prev.option <- options()$show.signif.stars
	options(show.signif.stars = FALSE)

	cat("ANOVA TABLE\nResponse Variable: ", trimStrings(strsplit(mymodel, split = "~")[[1]][1]), "\n\n", sep = "")
	printAOVTable(AOVContrastTable)
	cat("\n")
     
	contrast.option$coef <- cm
	options(show.signif.stars = prev.option)
	invisible(list(anova.contrast = AOVContrastTable, contrast.option = contrast.option))
}### end stmt --- contrast.compute
