# confInt2 - obtained from package nlrwr

confint2 <- function (object, parm, level = 0.95, method = c("asymptotic", 
                                                 "profile"), ...) 
{
  method <- match.arg(method)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- seq_along(pnames)
  if (is.numeric(parm)) 
    parm <- pnames[parm]
  asCI <- function(object, parm, level) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- stats:::format.perc(a, 3)
    fac <- qt(a, df.residual(object))
    parmInd <- match(parm, pnames)
    ci <- array(NA, dim = c(length(parmInd), 2), dimnames = list(parm, 
                                                                 pct))
    ses <- sqrt(diag(vcov(object)))[parmInd]
    ci[] <- cf[parmInd] + ses %o% fac
    ci
  }
  asProf <- function(object, parm, level) {
    message("Waiting for profiling to be done...")
    utils::flush.console()
    object <- profile(object, which = parm, alphamax = (1 - 
                                                          level)/4)
    confint(object, parm = parm, level = level, ...)
  }
  switch(method, asymptotic = asCI(object, parm, level), profile = asProf(object, 
                                                                          parm, level))
}
# <environment: namespace:nlrwr>