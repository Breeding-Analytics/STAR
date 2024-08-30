#######mshapiro.test
###Performs Shapiro-Wilk Multivariate Normality Test
###lifted from the mvnormtest package

mshapiro.test <- function (U) 
{
  U <- na.omit(U)
  if (!is.matrix(U)) 
    stop("U[] is not a matrix with number of columns (sample size) between 3 and 5000")
  n <- ncol(U)
  if (n < 3 || n > 5000) 
    stop("sample size must be between 3 and 5000")
  rng <- range(U)
  rng <- rng[2] - rng[1]
  if (rng == 0) 
    stop("all `U[]' are identical")
  Us <- apply(U, 1, mean)
  R <- U - Us
  M.1 <- solve(R %*% t(R), tol = 1e-18)
  Rmax <- diag(t(R) %*% M.1 %*% R)
  C <- M.1 %*% R[, which.max(Rmax)]
  Z <- t(C) %*% U
  return(shapiro.test(Z))
}