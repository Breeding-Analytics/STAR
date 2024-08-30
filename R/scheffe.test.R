# ----------------------------
# FROM agricolae package
# with some modification in presenting the output
# ----------------------------

# AAGulles added the parameter order
`scheffe.test` <- function (y, trt, DFerror, MSerror, Fc, alpha=0.05, group=TRUE,main = NULL, pwOrder = c("trmt", "means")) {
	pwOrder <- match.arg(pwOrder)	# added by AAGulles
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
    	clase<-c("aov","lm")
    	if("aov"%in%class(y) | "lm"%in%class(y)){
    		A<-y$model
    		DFerror<-df.residual(y)
    		MSerror<-deviance(y)/DFerror
    		Fc<-anova(y)[trt,4]
    		y<-A[,1]
    		ipch<-pmatch(trt,names(A))
    		if( is.na(ipch)) return(cat("Name: ",trt,"\n",names(A)[-1],"\n"))
    		name.t <-names(A)[ipch]
    		trt<-A[,ipch]
    		name.y <- names(A)[1]
    	}
	if (group && length(unique(trt)) > 26 ) { group <- FALSE }
    	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    	means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
    	sds <-   tapply.stat(junto[,1],junto[,2],stat="sd") #change
    	nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
    	means<-data.frame(means,std.err=sds[,2]/sqrt(nn[,2]),replication=nn[,2])
    	names(means)[1:2]<-c(name.t,name.y)
	# row.names(means)<-means[,1]
    	ntr<-nrow(means)
    	Tprob <- qf(1-alpha,ntr-1, DFerror)
    	nr<- 1/mean(1/nn[,2])
	# AAGulles suppress the printing of the following: 
    	# cat("\nStudy:", main)
    	# cat("\n\nScheffe Test for",name.y,"\n") 
    	# cat("\nMean Square Error  :",MSerror,"\n\n")
    	# cat(paste(name.t,",",sep="")," means\n\n")
    	# print(data.frame(row.names = means[,1], means[,-1]))
	# cat("\nalpha:",alpha,"; Df Error:",DFerror,"\n")
	# cat("Critical Value of F:", Tprob,"\n")
	scheffe <- sqrt(Tprob*(ntr-1)*2*MSerror/nr)

	# added by AAGulles
	if (!is.null(main)) { cat("\n", main,"\n", sep = "") }
    	cat("\nScheffe Test\n\n", sep = "") # added by AAGulles
    	maxWidth <- max(nchar(round(MSerror)), nchar(DFerror), nchar(round(Tprob,0)), nchar(round(scheffe,0)), nchar(round(nr,0))) + 7
	labelWidth <- 29
    	cat(formatC("Alpha", format = "s", width = labelWidth, flag = "-"), formatC(alpha, format = "f", digits = 2, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Error Degrees of Freedom", format = "s", width = labelWidth, flag = "-"), formatC(DFerror, format = "d", width =  maxWidth, flag = "#"), "\n",sep = "")  	
    	cat(formatC("Error Mean Square", format = "s", width = labelWidth, flag = "-"), formatC(MSerror, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Critical Value", format = "s", width = labelWidth, flag = "-"), formatC(Tprob, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	cat(formatC("Test Statistics", format = "s", width = labelWidth, flag = "-"), formatC(scheffe, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
	if (length(unique(nn[,2]))!=1) {
		# cat("\nHarmonic Mean of Cell Sizes ", nr )
		cat(formatC("Harmonic Mean of Cell Sizes", format = "s", width = labelWidth, flag = "-"), formatC(nr, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
	}
	cat("\n")

	if (group) {
		# suppress printing by AAGulles
		# cat("\nMinimum Significant Difference:",scheffe,"\n")
		# cat("\nMeans with the same letter are not significantly different.")
		# cat("\n\nGroups, Treatments and means\n")
		output <- order.group(means[,1], means[,2], means[,4], MSerror, Tprob,means[,3], parameter=1)
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
	} else {
		comb <-combn(ntr,2)
		nn<-ncol(comb)
		dif<-rep(0,nn)
		sig<-NULL
		LCL<-dif
		UCL<-dif
		pvalue<-rep(0,nn)
		for (k in 1:nn) {
			i<-comb[1,k]
			j<-comb[2,k]	
			if (means[i, 2] < means[j, 2]){
				comb[1, k]<-j
				comb[2, k]<-i
			}
			dif[k]<-abs(means[i,2]-means[j,2])
			sdtdif<-sqrt(MSerror * (1/means[i,4] + 1/means[j,4]))
			pvalue[k]<- round(1-pf(dif[k]^2/((ntr-1)*sdtdif^2),ntr-1,DFerror),6)

			LCL[k] <- dif[k] - sqrt(Tprob*(ntr-1))*sdtdif
			UCL[k] <- dif[k] + sqrt(Tprob*(ntr-1))*sdtdif
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
		output<-data.frame("Difference" = dif, pvalue=pvalue,sig,LCL,UCL)
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
		# print(output)
		# output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4],std.err=means[,3])
	}
	invisible(list(method = "Scheffe Test", tabValue = Tprob, testStat = scheffe, summary = output))
}

