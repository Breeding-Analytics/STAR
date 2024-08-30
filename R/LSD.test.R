# ----------------------------
# FROM agricolae package
# with some modification in presenting the output
# ----------------------------

# AAGulles added the parameter order
`LSD.test` <- function (y, trt, DFerror, MSerror, alpha = 0.05, p.adj = c("none",
    "holm", "hochberg", "bonferroni", "BH", "BY", "fdr"), group = TRUE,
    main = NULL, pwOrder = c("trmt", "means")) {

    	p.adj <- match.arg(p.adj)
    	pwOrder <- match.arg(pwOrder)	# added by AAGulles
    	clase<-c("aov","lm")
    	name.y <- paste(deparse(substitute(y)))
    	# name.t <- paste(deparse(substitute(trt)))
	if (is.data.frame(trt)) { 
		name.t <- names(trt) 
		
	} else { 
		name.t <- paste(deparse(substitute(trt))) 
		index <- which(strsplit(noquote(name.t), split = "")[[1]] == "\"")
		if (length(index) == 2) {
			name.t <- substr(noquote(name.t), start = index[1], stop = index[2])
		}
	} # modified by AAGulles

    	if("aov"%in%class(y) | "lm"%in%class(y)){
    		A<-y$model
    		DFerror<-df.residual(y)
    		MSerror<-deviance(y)/DFerror
    		y<-A[,1]
    		ipch<-pmatch(trt,names(A))
    		if( is.na(ipch)) return(cat("Name: ",trt,"\n",names(A)[-1],"\n"))
    		name.t <-names(A)[ipch]
    		trt<-A[,ipch]
    		name.y <- names(A)[1]
    	}

	if (group && length(unique(trt)) > 26 ) { group <- FALSE }
    	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
	#cat(str(junto), "\n")
    	means <- tapply.stat(junto[, 1], junto[, 2], stat="mean") #change
    	sds <- tapply.stat(junto[, 1], junto[, 2], stat="sd")     #change
    	nn <- tapply.stat(junto[, 1], junto[, 2], stat="length")  #change
    	std.err <- sds[, 2]/sqrt(nn[, 2])
    	Tprob <- qt(1 - alpha/2, DFerror)
    	LCL <- means[,2]-Tprob*std.err
    	UCL <- means[,2]+Tprob*std.err
    	means <- data.frame(means, std.err, replication = nn[, 2],LCL, UCL)
    	names(means)[1:2] <- c(name.t, name.y)
    	#row.names(means) <- means[, 1]
    	ntr <- nrow(means)
    	nk <- choose(ntr, 2)
    	if (p.adj != "none") {
        	a <- 1e-06
        	b <- 1
        	for (i in 1:100) {
            	x <- (b + a)/2
            	d <- p.adjust(x, n = nk, p.adj) - alpha
            	fa <- p.adjust(a, n = nk, p.adj) - alpha
            	if (d * fa < 0)
                		b <- x
            	if (d * fa > 0)
                		a <- x
        	}
        	Tprob <- qt(1 - x/2, DFerror)
    	}
    	nr <- unique(nn[, 2])
    	
	# AAGulles suppress the printing of the following: 
    	#cat("\nStudy:", main)
    	#cat("\n\nLSD t Test for", name.y, "\n")
    	#if (p.adj != "none")
   	#    cat("P value adjustment method:", p.adj, "\n")
    	#cat("\nMean Square Error: ",MSerror,"\n\n")
    	#cat(paste(name.t,",",sep="")," means and individual (",(1-alpha)*100,"%) CI\n\n")
    	#print(data.frame(row.names = means[,1], means[,-1]))
    	#cat("\nalpha:",alpha,"; Df Error:",DFerror)
    	#cat("\nCritical Value of t:", Tprob,"\n")
    
    	if (!is.null(main)) { cat("\n", main,"\n", sep = "") }
    	cat("\nLeast Significant Difference (LSD) Test\n\n", sep = "") # added by AAGulles
    	# the following was taken from the if(group)-else stmt and modified by AAGulles
    	if (length(nr) == 1) { 
		LSD <- Tprob * sqrt(2 * MSerror/nr)
		maxWidth <- max(nchar(round(MSerror)), nchar(DFerror), nchar(round(Tprob,0)), nchar(round(LSD,0))) + 7
    	} else {
        	nr1 <- 1/mean(1/nn[, 2])
        	LSD <- Tprob * sqrt(2 * MSerror/nr1)
		maxWidth <- max(nchar(round(MSerror)), nchar(DFerror), nchar(round(Tprob,0)), nchar(round(LSD,0)), nchar(round(nr1,0))) + 7
    	}

    	# AAGulles added the following:
    	labelWidth <- 29
    	if (p.adj != "none") { maxWidth <- max(nchar(p.adj), (maxWidth - 2)) + 2 }
    	cat(formatC("Alpha", format = "s", width = labelWidth, flag = "-"), formatC(alpha, format = "f", digits = 2, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Error Degrees of Freedom", format = "s", width = labelWidth, flag = "-"), formatC(DFerror, format = "d", width =  maxWidth, flag = "#"), "\n",sep = "")  	
    	cat(formatC("Error Mean Square", format = "s", width = labelWidth, flag = "-"), formatC(MSerror, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Critical Value", format = "s", width = labelWidth, flag = "-"), formatC(Tprob, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Test Statistics", format = "s", width = labelWidth, flag = "-"), formatC(LSD, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
	if (length(nr) > 1) {
		cat(formatC("Harmonic Mean of Cell Sizes", format = "s", width = labelWidth, flag = "-"), formatC(nr1, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
	}
    	if (p.adj != "none") {
  	   	cat(formatC("P value Adjustment Method", format = "s", width = labelWidth, flag = "-"), formatC(p.adj, format = "s", width = maxWidth, flag = "#"), "\n",sep = "")  
    	}
    	cat("\n")
	
    	if (group) {
	  	# suppress printing by AAGulles
        	#if (length(nr) == 1) {
        	#    LSD <- Tprob * sqrt(2 * MSerror/nr)
	  	#    cat("\nLeast Significant Difference", LSD) 
        	#} else {
        	#    nr1 <- 1/mean(1/nn[, 2])
	  	#	 LSD1 <- Tprob * sqrt(2 * MSerror/nr1)
        	#    cat("\nLeast Significant Difference", LSD1)
        	#    cat("\nHarmonic Mean of Cell Sizes ", nr1)
       	#}
        	#cat("\nMeans with the same letter are not significantly different.")
        	#cat("\n\nGroups, Treatments and means\n")
        	output <- order.group(means[, 1], means[, 2], means[, 4], MSerror, Tprob, means[, 3])
	  	colnames(output)[1] <- name.t  # added by AAGulles
	 	# AAGulles added the following if-else stmt:
	  	if (pwOrder == "means") { w<-order(means[,2],decreasing = TRUE) ## original LSD.test
	  	} else {
	  	     if (suppressWarnings(all(!is.na(as.numeric(as.character(factor(trimStrings(output[,1])))))))) {
	  	          output[,1] <- factor(as.numeric(as.character(output[,1])))
	  	     }
			output <- output[order(output[,1]),]
			rownames(output) <- 1:nrow(output)
			w<-order(means[,1])
        	}
	  	cat("Summary of the Result:\n", sep = "") # added by AAGulles
	  	printDataFrame(output[,c(1,2,3,5)])	  # added by AAGulles
	  	cat("Means with the same letter are not significantly different.\n\n") # added by AAGulles
        	output <- data.frame(output,LCI=means[w,5],UCI=means[w,6])
    	}

    	if (!group) {
        	comb <- combn(ntr, 2)
        	nn <- ncol(comb)
        	dif <- rep(0, nn)
        	LCL1<-dif
		UCL1<-dif
        	sig<-NULL
        	pvalue <- rep(0, nn)
        	for (k in 1:nn) {
            	i <- comb[1, k]
            	j <- comb[2, k]
            	if (means[i, 2] < means[j, 2]){
            		comb[1, k]<-j
            		comb[2, k]<-i
            	}
            	dif[k] <- abs(means[i, 2] - means[j, 2])
           	 	sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,4]))
            	pvalue[k] <- 2 * (1 - pt(dif[k]/sdtdif, DFerror))
            	if (p.adj != "none") pvalue[k] <- p.adjust(pvalue[k], n = nk, p.adj)
            	pvalue[k] <- round(pvalue[k],6)
        		LCL1[k] <- dif[k] - Tprob*sdtdif
			UCL1[k] <- dif[k] + Tprob*sdtdif
        		sig[k]<-" "
	  		# AAGulles suppress the following if-else stmt:
        		#if (pvalue[k] <= 0.001) sig[k]<-"***"
        		#else  if (pvalue[k] <= 0.01) sig[k]<-"**"
        		#else  if (pvalue[k] <= 0.05) sig[k]<-"*"
        		#else  if (pvalue[k] <= 0.1) sig[k]<-"."

			# added by AAGulles
			remark <- "  "					
			if (alpha <= 0.001)  { remark <-"***"	
        		} else {
				if (alpha <= 0.01) { remark <-"**"	
				} else {
					if (alpha <= 0.05) { remark<-"*"	
					} else {
						if (alpha <= 0.1) remark <-"."	
					}
				}
			}  
			if (pvalue[k] <= alpha) { sig[k] <- remark } else { sig[k] <- "   " }
        	}
        	tr.i <- means[comb[1, ],1]
        	tr.j <- means[comb[2, ],1]
        	output<-data.frame("Difference" = dif, pvalue = pvalue,sig,LCL=LCL1,UCL=UCL1)
        	rownames(output)<-paste(tr.i,tr.j,sep=" - ")
		output <- output[,c("Difference", "pvalue", "sig")]
		
        	#cat("Comparison between treatments means\n")
	  	#printDataFrame(cbind("Treatment" = rownames(output),output))	# change by AAGulles
		sigResult <- output[output[,"sig"] == remark, c("Difference", "pvalue")]
		colnames(output) <- c("MeanDiff", "Prob", "Sig")
		if (nrow(sigResult) != 0) {
			colnames(sigResult) <- c("Mean Diff", "Prob")
			cat("Significant Pairwise Mean Comparison at alpha = ", alpha,"\n")
		  	printDataFrame(cbind("Treatment" = rownames(sigResult),sigResult))	# change by AAGulles
		}
	 
	 	# AAGulles suppress the following:
        	# output <- data.frame(trt = means[, 1], means = means[,2], M = "", N = means[, 4], std.err ,LCL,UCL)
    	}
    	invisible(list(method = "Least Significant Difference (LSD) Test", tabValue = Tprob, testStat = LSD, summary = output))
}
