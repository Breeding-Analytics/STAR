#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#	CLUSTER ANALYSIS FUNCTION (AGGLOMERATIVE) 09.25.2012 Modified: 11.08.13                      								     			     			                                                                                       
#	ClusterAgglo <- function( data,var,idVar,sbinVar,abinVar,ofactorVar,factorVar,stand= TRUE,distance ,clusmethod,distMatrix,
#                              clusterMem,descriptiveStatRaw,corMatx,scatterMatx,descriptiveStat,dendrogram,clusterBox,clusterNum,
#                              saveMem,outputPath)
#		       
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

ClusterAgglo <- function(data, var = NULL, idVar = NULL, sbinVar = NULL, abinVar = NULL, ofactorVar = NULL,factorVar = NULL, stand= TRUE, 
                         distance = c("Euclidean", "Maximum", "Manhattan", "Minkowski", "Canberra", "Binary", "Simple Matching","Sokal & Sneath", "Hamann coefficient","Jaccard", "Dice", "Gower"), 
                         clusmethod = c("Single", "Complete", "Average", "Ward", "Centroid"),distMatrix = FALSE, copDistance = FALSE, clusterMem = TRUE, 
                         descriptiveStatRaw = TRUE, corMatx = TRUE, scatterMatx = TRUE, descriptiveStat = TRUE, dendrogram = TRUE, clusterBox = TRUE, clusterNum=2, saveMem = TRUE, outputPath = NULL) 
  UseMethod("ClusterAgglo")  

ClusterAgglo.default <- function(data, var = NULL, idVar = NULL, sbinVar = NULL, abinVar = NULL, ofactorVar = NULL,factorVar = NULL, stand= TRUE, 
                                 distance = c("Euclidean", "Maximum", "Manhattan", "Minkowski", "Canberra", "Binary", "Simple Matching","Sokal & Sneath", "Hamann coefficient","Jaccard", "Dice", "Gower"), 
                                 clusmethod = c("Single", "Complete", "Average", "Ward", "Centroid"),distMatrix = FALSE,copDistance = FALSE, clusterMem = TRUE, 
                                 descriptiveStatRaw = TRUE, corMatx = TRUE, scatterMatx = TRUE, descriptiveStat = TRUE, dendrogram = TRUE, clusterBox = TRUE, clusterNum=2, saveMem = TRUE, outputPath = NULL) {
  if (is.character(data)) { 
    nameData <- data
    if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- eval(parse(text = data)) 
  } else {
    if (is.data.frame(data)) { 
      nameData <- paste(deparse(substitute(data)))  
      tempData <- data  
    } else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
  } 
  
  numObsNM = nrow((na.omit(tempData[,c(var,idVar,sbinVar,abinVar,ofactorVar,factorVar)])))
  numObs = nrow((tempData[,c(var,idVar,sbinVar,abinVar,ofactorVar,factorVar)]))
  
  if (numObs == numObsNM) { cat("Number of Observations:", numObs,"\n")
  } else {
    cat("Number of Observations:", numObs,"\n")
    cat("Number of Observations Used:", numObsNM,"\n\n")
  }
  
  tempData <- na.omit(tempData[,c(var,idVar,sbinVar,abinVar,ofactorVar,factorVar)])
  
  if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
  
  if (!is.null(var))     { 
    if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }  
    if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
    tempVar <- tempData[var]  
  }
  
  if (!is.null(idVar))   	{ 
    if (any(is.na(match(idVar, names(tempData))))) { stop("At least one item in the character vector 'idVAr' does not match any variable name in the dataset.") }
    tempGrp <- tempData[idVar]
    for (i in (1:ncol(tempGrp))) { tempGrp[,i] <- factor(tempGrp[,i]) }
  } else { tempGrp <- rep(1, each = nrow(tempData)) }
  
  if (!is.null(sbinVar))     { 
    if (any(is.na(match(sbinVar, names(tempData))))) { stop("At least one item in the character vector 'sbinVar' does not match any variable name in the dataset.") }
    tempsBinVar <- tempData[sbinVar]
    for (i in (1:ncol(tempsBinVar))) { tempsBinVar[,i] <- factor(tempsBinVar[,i]) }
    if(nlevels(tempsBinVar[,i] != 2)) { stop("The variable is not binary") }
  } else { tempsBinVar <- rep(1, each = nrow(tempData)) }
  
  if (!is.null(abinVar))     { 
    if (any(is.na(match(abinVar, names(tempData))))) { stop("At least one item in the character vector 'abinVar' does not match any variable name in the dataset.") }
    tempaBinVar <- tempData[abinVar]
    for (i in (1:ncol(tempaBinVar))) { tempaBinVar[,i] <- factor(tempaBinVar[,i]) }
    if(nlevels(tempaBinVar[,i] != 2)) { stop("The variable is not binary") }
  } else { tempaBinVar <- rep(1, each = nrow(tempData)) }
  
  if (!is.null(ofactorVar))     { 
    if (any(is.na(match(ofactorVar, names(tempData))))) { stop("At least one item in the character vector 'ofactorVar' does not match any variable name in the dataset.") }
    tempofactorVar <- tempData[ofactorVar]
    for (i in (1:ncol(tempofactorVar))) { tempofactorVar[,i] <- as.ordered(tempofactorVar[,i]) }
  } else { tempofactorVar <- rep(1, each = nrow(tempData)) }
  
  
  if (!is.null(factorVar))     { 
    if (any(is.na(match(factorVar, names(tempData))))) { stop("At least one item in the character vector 'factorVar' does not match any variable name in the dataset.") }
    tempfactorVar <- tempData[factorVar]
    for (i in (1:ncol(tempfactorVar))) { tempfactorVar[,i] <- factor(tempfactorVar[,i]) }
  } else { tempfactorVar <- rep(1, each = nrow(tempData)) }
  
  distance <- match.arg(distance)
  clusmethod <- match.arg(clusmethod)
  options(width = 6000, digits =3)
  
  # ------------------first condition var not null---------------------- 
  
  if(!is.null(var) & (is.null(sbinVar) & is.null(abinVar) & is.null(ofactorVar) & is.null(factorVar))){
    x <- tempData[,var]
    z <- var
    a <- data.Normalization(x, type = "n1")
    if(stand){
      adist <- dist(a, method = tolower(distance))
    } else adist <- dist(x, method = tolower(distance))
  }
  
  #---------------------------------------------------------------------  
  #------------------second condition binVar not null-------------------  
  
  if(!is.null(sbinVar) & (is.null(var) & is.null(abinVar) & is.null(ofactorVar) & is.null(factorVar))){
    x<- tempData[,sbinVar]
    z <- sbinVar
    if(distance == "Simple Matching"){
      adist <- dist.binary(x, method = 2)
    }
    if(distance == "Sokal & Sneath"){
      adist <- dist.binary(x, method = 3)
    }
    if(distance == "Hamann coefficient"){
      adist <- dist.binary(x, method = 6)
    }
  }
  
  if(!is.null(abinVar) & (is.null(var) & is.null(sbinVar) & is.null(ofactorVar) & is.null(factorVar))){
    x<- tempData[,abinVar]
    z <- abinVar
    if(distance == "Jaccard"){
      adist <- dist(x, method = "binary")
    }
    if(distance == "Dice"){
      adist <- dist.binary(x, method = 5)
    }   
  }
  
  if(!is.null(factorVar) & (is.null(var) & is.null(sbinVar) & is.null(ofactorVar) & is.null(abinVar))){
    x<- tempData[,factorVar]
    z<- factorVar
    if(stand){   
      adist <- daisy(x, stand=TRUE, metric = tolower(distance))
    }else adist <- daisy(x, metric = tolower(distance))
  } 
  
  #--------------------------------------------------------------------  
  #------------------third condition all not null----------------------
  
  if((!is.null(var) & !is.null(sbinVar)) & is.null(abinVar) & is.null(factorVar) & is.null(ofactorVar)){
    x<- tempData[c(var,sbinVar)]
    z<- c(var,sbinVar)
    if(stand){   
      #     adist <- daisy(x, stand = TRUE, type = list(symm = names(tempData[sbinVar])))
      #       adist <- daisy(x, stand = TRUE, type = list(symm = tempData[,sbinVar]))
      adist <- daisy(x, stand = TRUE, type = list(symm = sbinVar))
    }else 
      #     adist <- daisy(x, type = list(symm = names(tempData[sbinVar])))
      #       adist <- daisy(x,type = list(symm = tempData[,sbinVar]))
      adist <- daisy(x, type = list(symm = sbinVar))
  }
  
  if((!is.null(var) & !is.null(abinVar)) &  is.null(sbinVar)& is.null(factorVar) & is.null(ofactorVar)){
    x<- tempData[c(var,abinVar)]
    z<- c(var,abinVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(asymm = abinVar))
    }else adist <- daisy(x, type = list(asymm = abinVar))
  }
  if((!is.null(var) & !is.null(ofactorVar)) &  is.null(sbinVar)& is.null(factorVar) & is.null(abinVar)){
    x<- tempData[c(var,ofactorVar)]
    z<- c(var,ofactorVar)
    if(stand){
      adist <- daisy(x, stand = TRUE, type = list(ordratio = ofactorVar))
    }else  adist <- daisy(x, type = list(ordratio = ofactorVar))
  }
  if((!is.null(var) & !is.null(factorVar)) & is.null(sbinVar) & is.null(abinVar) & is.null(ofactorVar)){
    x<- tempData[c(var,factorVar)]
    z<- c(var, factorVar) 
    if(stand){
      adist <- daisy(x,stand = TRUE, metric = tolower(distance))
    } else  adist <- daisy(x, metric = tolower(distance))
  }      
  
  
  
  if((!is.null(sbinVar) & !is.null(abinVar)) &  is.null(var)& is.null(factorVar) & is.null(ofactorVar)){
    x<- tempData[c(sbinVar,abinVar)]
    z<- c(sbinVar,abinVar)
    if(stand){ 
      adist <- daisy(x, stand = TRUE, type = list(asymm = abinVar, symm = sbinVar))
      #     adist <- daisy(x, stand = TRUE, type = list(asymm = as.vector(names(tempData[,abinVar])), symm = as.vector(names(tempData[,sbinVar]))))
    }else 
      #       adist <- daisy(x, type = list(asymm = as.vector(names(tempData[,abinVar])), symm = as.vector(names(tempData[,sbinVar]))))
      adist <- daisy(x, type = list(asymm = abinVar, symm = sbinVar))
  }
  if((!is.null(sbinVar) & !is.null(ofactorVar)) &  is.null(var)& is.null(factorVar) & is.null(abinVar)){
    x<- tempData[c(sbinVar,ofactorVar)]
    z<- c(sbinVar,ofactorVar)
    if(stand){    
      adist <- daisy(x, stand= TRUE, type = list(ordratio = ofactorVar, symm = sbinVar))
    }else adist <- daisy(x, type = list(ordratio = ofactorVar, symm = sbinVar))
  }
  if((!is.null(sbinVar) & !is.null(factorVar)) &  is.null(ofactorVar)& is.null(var) & is.null(abinVar)){
    x<- tempData[c(sbinVar,factorVar)]
    z<- c(sbinVar,factorVar)
    if(stand){  
      adist <- daisy(x,stand=TRUE, type = list(symm = sbinVar))
    } else adist <- daisy(x, type = list(symm = sbinVar))
  }
  
  
  if((!is.null(abinVar) & !is.null(ofactorVar)) &  is.null(var)& is.null(factorVar) & is.null(sbinVar)){
    x<- tempData[c(abinVar,ofactorVar)]
    z<- c(abinVar,ofactorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(ordratio = ofactorVar, asymm = abinVar))
    }else adist <- daisy(x, type = list(ordratio = ofactorVar, asymm = abinVar))
  }
  if((!is.null(abinVar) & !is.null(factorVar)) &  is.null(ofactorVar)& is.null(var) & is.null(sbinVar)){
    x<- tempData[c(abinVar,factorVar)]
    z<- c(abinVar,factorVar)
    if(stand){
      adist <- daisy(x, stand = TRUE,type = list(asymm = abinVar))
    }else adist <- daisy(x, type = list(asymm = abinVar))
  }
  
  if((!is.null(factorVar) & !is.null(ofactorVar)) &  is.null(var)& is.null(sbinVar) & is.null(abinVar)){
    x<- tempData[c(factorVar,ofactorVar)]
    z<- c(factorVar,ofactorVar)
    if(stand){ 
      adist <- daisy(x,stand=TRUE, type = list(ordratio = ofactorVar))
    }else   adist <- daisy(x, type = list(ordratio = ofactorVar))
  }
  
  # three entries not null
  
  if((!is.null(var) & !is.null(sbinVar) & !is.null(abinVar)) & is.null(ofactorVar) & is.null(factorVar)){
    x<- tempData[c(var,sbinVar,abinVar)]
    z<- c(var,sbinVar,abinVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(symm = sbinVar, asymm = abinVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar))
  }
  if((!is.null(var) & !is.null(sbinVar) & !is.null(ofactorVar)) & is.null(abinVar) & is.null(factorVar)){
    x<- tempData[c(var,sbinVar,ofactorVar)]
    z<- c(var,sbinVar,ofactorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(symm = sbinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, ordratio = ofactorVar))
  }
  if((!is.null(var) & !is.null(sbinVar) & !is.null(factorVar)) & is.null(ofactorVar) & is.null(abinVar)){
    x<- tempData[c(var,sbinVar,factorVar)]
    z<- c(var,sbinVar,factorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(symm = sbinVar))
    }else adist <- daisy(x, type = list(symm = sbinVar))
  }
  if((!is.null(var) & !is.null(abinVar) & !is.null(ofactorVar)) & is.null(factorVar) & is.null(sbinVar)){
    x<- tempData[c(var,abinVar,ofactorVar)]
    z<- c(var,abinVar,ofactorVar)
    if(stand){   
      adist <- daisy(x, stand = TRUE,type = list(asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(asymm = abinVar, ordratio = ofactorVar))
  }
  if((!is.null(var) & !is.null(abinVar) & !is.null(factorVar)) & is.null(ofactorVar) & is.null(sbinVar)){
    x<- tempData[c(var,abinVar,factorVar)]
    z<- c(var,abinVar,factorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(asymm= abinVar))
    }else adist <- daisy(x, type = list(asymm= abinVar))
  }
  if((!is.null(var) & !is.null(ofactorVar) & !is.null(factorVar)) & is.null(abinVar) & is.null(sbinVar)){
    x<- tempData[c(var,ofactorVar,factorVar)]
    z<- c(var,ofactorVar,factorVar)
    if(stand){    
      adist <- daisy(x,stand = TRUE, type = list(ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(ordratio = ofactorVar))
  }
  if((!is.null(sbinVar) & !is.null(abinVar) & !is.null(ofactorVar)) & is.null(var) & is.null(factorVar)){
    x<- tempData[c(sbinVar,abinVar,ofactorVar)]
    z<- c(sbinVar,abinVar,ofactorVar)
    if(stand){ 
      adist <- daisy(x,stand = TRUE, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
    }else  adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
  }    
  if((!is.null(sbinVar) & !is.null(abinVar) & !is.null(factorVar)) & is.null(var) & is.null(ofactorVar)){
    x<- tempData[c(sbinVar,abinVar,factorVar)]
    z<- c(sbinVar,abinVar,factorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE, type = list(symm = sbinVar, asymm = abinVar))
    }else  adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar))
  }
  if((!is.null(sbinVar) & !is.null(ofactorVar) & !is.null(factorVar)) & is.null(var) & is.null(abinVar)){
    x<- tempData[c(sbinVar,ofactorVar,factorVar)]
    z<- c(sbinVar,ofactorVar,factorVar)
    if(stand){
      adist <- daisy(x,stand = TRUE,type = list(symm = sbinVar, ordratio = ofactorVar))
    }else  adist <- daisy(x, type = list(symm = sbinVar, ordratio = ofactorVar))
  }
  if((!is.null(abinVar) & !is.null(ofactorVar) & !is.null(factorVar)) & is.null(var) & is.null(sbinVar)){
    x<- tempData[c(abinVar,ofactorVar,factorVar)]
    z<- c(abinVar,ofactorVar,factorVar)
    if(stand){   
      adist <- daisy(x, stand = TRUE, type = list(asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(asymm = abinVar, ordratio = ofactorVar))
  }
  
  #four entries not null
  if((!is.null(var) & !is.null(sbinVar) & !is.null(abinVar) & !is.null(ofactorVar))  & is.null(factorVar)){
    x<- tempData[c(var,sbinVar,abinVar,ofactorVar)]
    z<- c(var,sbinVar,abinVar,ofactorVar)
    if(stand){   
      adist <- daisy(x,stand= TRUE,type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
  }
  if((!is.null(var) & !is.null(sbinVar) & !is.null(factorVar) & !is.null(abinVar))  & is.null(ofactorVar)){
    x<- tempData[c(var,sbinVar,factorVar,abinVar)]
    z<- c(var,sbinVar,factorVar,abinVar)
    if(stand){   
      adist <- daisy(x,stand= TRUE, type = list(symm = sbinVar, asymm = abinVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar))
  }
  if((!is.null(var) & !is.null(sbinVar) & !is.null(factorVar) & !is.null(ofactorVar))  & is.null(abinVar)){
    x<- tempData[c(var,sbinVar,factorVar,ofactorVar)]
    z<- c(var,sbinVar,factorVar,ofactorVar)
    if(stand){    
      adist <- daisy(x,stand= TRUE, type = list(symm = sbinVar,ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(symm = sbinVar,ordratio = ofactorVar))
  } 
  if((!is.null(factorVar) & !is.null(abinVar) & !is.null(ofactorVar) & !is.null(var))  & is.null(sbinVar)){
    x<- tempData[c(factorVar,abinVar,ofactorVar,var)]
    z<- c(factorVar,abinVar,ofactorVar,var)
    if(stand){  
      adist <- daisy(x,stand= TRUE, type = list(asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(asymm = abinVar, ordratio = ofactorVar))
  }    
  if((!is.null(sbinVar) & !is.null(abinVar) & !is.null(factorVar) & !is.null(ofactorVar))  & is.null(var)){
    x<- tempData[c(sbinVar,abinVar,factorVar,ofactorVar)]
    z<- c(sbinVar,abinVar,factorVar,ofactorVar)
    if(stand){  
      adist <- daisy(x,stand= TRUE, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
  }
  
  #five entries not null
  if(!is.null(sbinVar) & !is.null(abinVar) & !is.null(factorVar) & !is.null(ofactorVar)  & !is.null(var)){
    x<- tempData[c(sbinVar,abinVar,factorVar,ofactorVar,var)]
    z<- c(var,sbinVar,abinVar,ofactorVar,factorVar)
    if(stand){ 
      adist <- daisy(x,stand= TRUE, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
    }else adist <- daisy(x, type = list(symm = sbinVar, asymm = abinVar, ordratio = ofactorVar))
  }
  
  #output of the analysis
  cat("\nAGGLOMERATIVE CLUSTER ANALYSIS\n")
  cat("\nSPECIFICATIONS:\n")
  cat("\t\tDistance Method:\t\t",distance,"\n",sep = "")
  cat("\t\tClustering Method:\t\t", clusmethod,"\n",sep = "")
  cat("\t\tNumber of Clusters:\t\t", clusterNum,"\n\n",sep = "")
  
  
  attr(adist, "Labels") <- tempData[,idVar]
  
  if (distMatrix){ write.csv(as.matrix(adist), row.names = TRUE, file=paste(outputPath, "DistanceMatrix.csv", sep = "")) }
  
  amethod <- hclust(adist, method = tolower(clusmethod)) 
  Membership <- cutree(amethod, k = as.numeric(clusterNum))
  names(Membership) <- tempData[,idVar]
  memData <- data.frame(Membership)
  memberList <-list()
  
  if (descriptiveStatRaw){
    cat("\n")
    DescriptiveStatistics(data = tempData, var = var, grp = NULL, statistics = c("mean", "sd", "min", "max") )
    
  }
  
  if (corMatx){
    cat("\nCORRELATION MATRIX\n")
    cat("\n")
    cor(tempVar, method = "pearson")
    outCor <- round(cor(tempVar, method = "pearson"),4)
    print(outCor)
    cat("\n")
  }
  
  if (descriptiveStat){
    cat("\n")
    cat("\n")
    Cluster <- Membership
    all <- cbind(x, Cluster )
    DescriptiveStatistics(data = all, var = var, grp = "Cluster", statistics = c("mean", "sd", "min", "max") )
  } 
  
  
  if (clusterMem){  
    cat("\nCLUSTER MEMBERSHIP SUMMARY\n")
    memSummary <- cbind(tempData[idVar], x, memData)
    
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
    
    cat("\nNumber of members in each cluster\n")
    a<- data.frame(table(Membership))
    colnames(a)<- "Cluster"
    colnames(a)[2] <-"Size"
    printDataFrame(a)
  }
  #display plot and dendrogram
  if (!is.null(outputPath)) {
    if(scatterMatx){
      png(filename = paste(outputPath,"ScatterPlotMatrix.png",sep=""))  
      pairs(tempData[,c(var,sbinVar,abinVar,ofactorVar,factorVar)],main="ScatterPlot Matrix")
      dev.off()
    }
    if (dendrogram){
      png(paste(outputPath,"AggloGraph.png", sep=""),width=2000,height=1000)
      par(cex=0.8,font=3)
      if(!is.null(idVar)){amethod$labels <- tempData[,idVar]}
      dendro <- as.dendrogram(amethod)
      plot(amethod)
      plot(dendro,main = "Dendrogram using Agglomerative Clustering Method")
      
      if(clusterBox){
        rect.hclust(tree=amethod, k= as.numeric(clusterNum), border=c("red", "blue", "green", "purple"))
      }
      dev.off()
    }
    if (saveMem){write.csv(data.frame(Rownames= rownames(memSummary), memSummary), row.names = FALSE, file=paste(outputPath, "MembershipSummary.csv", sep = ""))}
  }
  
  Copcorr <- cophenetic(amethod)
  cop <- cor(adist, Copcorr)   
  cat("\nCOPHENETIC CORRELATION COEFFICIENT = ", cop, "\t\n")
  
  if(copDistance){write.csv(as.matrix(cophenetic(amethod)), row.names = TRUE, file=paste(outputPath, "CopheneticDistances.csv", sep = "")) }
  
  return(list(ClusterMethod = amethod, Membership = memSummary, DistMatrix = adist))
}
#### end statement ClusterAgglo####
