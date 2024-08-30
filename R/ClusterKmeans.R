#------------------------------------------------------------------------------------------------------------------------------------#
#     CLUSTER ANALYSIS FUNCTION (KMEANS)09.25.2012 MODIFIED: 03.15.13     				                       
#     ClusterKmeans <- function(data, var, clusterMem, descriptiveStatRaw, corMatx, descriptiveStat, kgraph, clusterNum, outputPath)  
#------------------------------------------------------------------------------------------------------------------------------------#

ClusterKmeans <- function(data, var, idVar= NULL, clusterMem = TRUE, descriptiveStatRaw = TRUE, corMatx = TRUE, descriptiveStat= TRUE,kgraph = TRUE, clusterNum=2, saveMem = TRUE, outputPath = NULL)
  
  UseMethod("ClusterKmeans")

ClusterKmeans.default <- function(data, var, idVar= NULL, clusterMem = TRUE, descriptiveStatRaw = TRUE, corMatx = TRUE, descriptiveStat= TRUE,kgraph = TRUE, clusterNum=2, saveMem = TRUE, outputPath = NULL)
{
  if (is.character(data)) { 
    nameData <- data
    if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- eval(parse(text = data)) 
  }else{
    if (is.data.frame(data)) { 
      nameData <- paste(deparse(substitute(data)))	
      tempData <- data	
    } else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
  }
  numObsNM = nrow((na.omit(tempData[,c(var,idVar)])))
  numObs = nrow((tempData[,c(var,idVar)]))
  
  if (numObs == numObsNM) { cat("Number of Observations:", numObs,"\n")
  } else {
    cat("Number of Observations:", numObs,"\n")
    cat("Number of Observations Used:", numObsNM,"\n\n")
  }
  
  tempData <- na.omit(tempData[,c(var,idVar)])
  
  if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
  if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }
  if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
  
  if (!is.null(idVar)) { 
    if (!is.character(idVar)) { stop(paste("The object 'idVar' should be a character vector.", sep = "")) }
    if (any(is.na(match(idVar, names(tempData))))) { stop("At least one item in the character vector 'idVar' does not match any variable name in the dataset.") }
    for (i in (1:length(idVar))) { if (!is.factor(tempData[,idVar[i]])) { tempData[,idVar[i]] <- factor(tempData[,idVar[i]]) }}
  }
  
  options(width = 6000, digits = 3)
  
  #output of the analysis
  km <- kmeans(tempData[,var], centers = as.numeric(clusterNum))
  
  cat("K-MEANS CLUSTER ANALYSIS\n")
  
  Membership <- km$cluster
  names(Membership) <- tempData[,idVar]
  memData <- data.frame(Membership)
  memberList <- list()
  
  if (descriptiveStatRaw){
    cat("\n")
    DescriptiveStatistics(data = tempData, var = var, grp = NULL, statistics = c("mean", "sd", "min", "max") )
    
  }
  
  if (corMatx){
    cat("\nCORRELATION MATRIX\n")
    cat("\n")
    cor(tempData[,var], method = "pearson")
    outCor <- round(cor(tempData[,var], method = "pearson"),4)
    print(outCor)
    cat("\n")
  }
  
  if (descriptiveStat){
    Cluster <- Membership
    all <- cbind(tempData[,var], Cluster)
    DescriptiveStatistics(data = all, var = var, grp = "Cluster", statistics = c("mean", "sd", "min", "max") )
  }
  
  
  if (clusterMem){
    cat("\nCLUSTER MEMBERSHIP SUMMARY\n")
    cat("\t\tNumber of Clusters:\t\t", clusterNum,"\n",sep = "" )	
    memSummary <- cbind(tempData[idVar], tempData[,var], memData)
    for (i in (1:as.numeric(clusterNum))) {
      cat("\nMember of Cluster ",i,"\n")
      temp <- rownames(subset(memData, Membership == i))
      memberList[[i]] <- rownames(subset(memData, Membership == i))
      names(memberList)[i] <- paste("Cluster Number:", i)
      index <- 1
      for (j in (1:ceiling(length(temp)/15))) {
        if(index+14 > length(temp)) { cat(temp[index:length(temp)], "\n")
        } else { cat(temp[index:(index+14)], "\n") }
        index <- index + 15
      }
    }
    cat("\n")
    
    
    cat("\nNumber of members in each cluster\n")
    a<- data.frame(table(Membership))
    colnames(a)<- "Cluster"
    colnames(a)[2] <-"Size"
    printDataFrame(a)
    
    # 	    	cat("\nNumber of members in each cluster\n\n")
    #     		cat("Cluster number\t")
    # 	    	for (i in (1:attributes(table(Membership))$dim)) {
    #     		  	cat("\t",attributes(table(Membership))$dimnames$Membership[i])
    # 	    	}
    #     		cat("\n")
    # 	    	cat("Cluster Size\t")
    #     		for (i in (1:attributes(table(Membership))$dim)) {
    #       		cat("\t",table(Membership)[i])
    # 	    	}
  }
  cat("\n")
  cat("\n")
  
  
  #display plot
  if (!is.null(outputPath)) {
    if (kgraph){
      png(paste(outputPath, "Kmeans.png", sep = ""))
      plot(tempData[,var], col=Membership, pch=Membership)
      dev.off()
    }
    if (saveMem){
      write.csv(data.frame(Rownames= rownames(memSummary), memSummary), row.names = FALSE, file=paste(outputPath, "MembershipSummary.csv", sep = ""))
    }
  }
  return(list(ClusterMethod = km, Membership = memSummary))
  
}#### end statement ClusterKmeans####