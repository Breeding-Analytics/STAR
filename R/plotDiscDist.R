# --------------------------------------------------------       
# plotDiscDist: plots empirical and theoretical discrete distributions
#               modified "plotdist" function of the fitdistrplus package
#
# ARGUMENTS:
# data - name of R dataframe
# distr - discrete distribution(s) to fit
# para - list of parameters for the fitted distribution(s)
# --------------------------------------------------------

plotDiscDist <- function (data, distr, para, ...) 
  UseMethod("plotDiscDist")

plotDiscDist.default <- function (data, distr, para, ...) 
{
  def.par <- par(no.readonly = TRUE)
  if (missing(data) || !is.vector(data, mode = "numeric")) 
    stop("data must be a numeric vector")
  if ((missing(distr) & !missing(para)) || (missing(distr) & 
                                              !missing(para))) 
    stop("distr and para must defined")
  xlim <- c(min(data), max(data))

  if (!is.character(distr)) {
    distname <- substring(as.character(match.call()$distr), 2)
  } else distname <- distr
  if (!is.list(para)) 
    stop("'para' must be a named list")
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  pdistname <- paste("p", distname, sep = "")
  if (!exists(pdistname, mode = "function")) 
    stop(paste("The ", pdistname, " function must be defined"))
  qdistname <- paste("q", distname, sep = "")
  if (!exists(qdistname, mode = "function")) 
    stop(paste("The ", qdistname, " function must be defined"))
  densfun <- get(ddistname, mode = "function")
  nm <- names(para)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m))) 
    stop(paste("'para' specifies names which are not arguments to ", ddistname))
  n <- length(data)
  par(mfrow = c(1, 2))
  t <- table(data)
  xval <- as.numeric(names(t))
  xvalfin <- seq(min(xval), max(xval))
  xlinesdec <- min((max(xval) - min(xval))/30, 0.4)

  numDistr = length(distr)
  ydobs <- as.vector(t)/n
  yd <- NULL
  ydmax <- max(ydobs)
  plotColors <- c("red", "blue", "green")
  lineTypes <- NULL
  for (i in 1:numDistr) {
    yd[[i]] <- do.call(ddistname[[i]], c(list(x = xvalfin), as.list(para[[i]]$estimate)))
    ydmax <- max(ydmax, yd[[i]])
    lineTypes[[i]] = i + 1
  }  
  plot(xvalfin + xlinesdec, yd[[1]], type = "h", xlim = c(min(xval), 
       max(xval) + xlinesdec), ylim = c(0, ydmax), lty = lineTypes[1], 
       col = plotColors[1], main = "Emp. and theo. distr.", 
       xlab = "Data", ylab = "Density", ...)
  if (numDistr > 1)
    for (i in 2:numDistr)
      points(xvalfin + xlinesdec, yd[[i]], type = "h", lty = lineTypes[i], col = plotColors[i], ...)
  lineTypes[[numDistr+1]] <- 2
  plotColors[[numDistr+1]] <- "black"
  legLabel <- c(distr, "empirical")
  points(xval, ydobs, type = "h", lty = lineTypes[[numDistr+1]], col = plotColors[[numDistr+1]], ...)
  legend("topright", lty = lineTypes, col = plotColors, legend = legLabel, 
      bty = "o", bg = "white", cex = 0.6, ...)
  
  ycdfobs <- ecdf(data)(xvalfin)
  ycdf <- NULL

  for (i in 1:numDistr)
    ycdf[[i]] <- do.call(pdistname[[i]], c(list(q = xvalfin), as.list(para[[i]]$estimate)))
  
  plot(xvalfin + xlinesdec, ycdf[[1]], type = "h", xlim = c(min(xval), 
      max(xval) + xlinesdec), ylim = c(0, 1), lty = lineTypes[1], 
      col = plotColors[1], main = "Emp. and theo. CDFs", xlab = "Data", 
      ylab = "CDF", ...)
  if (numDistr > 1)  
    for (i in 2:numDistr)
      points(xvalfin + xlinesdec, ycdf[[i]], type = "h", lty = lineTypes[i], col = plotColors[i],  ...)
  points(xvalfin, ycdfobs, type = "h", lty = lineTypes[[numDistr+1]], col = plotColors[numDistr+1],  ...)
  legend("topleft", lty = lineTypes, col = plotColors, legend = legLabel, 
         bty = "o", bg = "white", cex = 0.6, ...)
  par(def.par)
}