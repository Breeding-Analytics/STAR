#----------------------------------------------------#
# BoxsMTest - function for performing Box's M Test
# Parameters: 
# X - data frame containing factor (grouping) and numeric variables
# vars - numeric variables
# factorVar - factor
# alpha - significance level (default = 0.05). 
#
# MODIFIED R FUNCTION OF ANDY LIAW 
#http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33330.html
# [R] evaluation of discriminant functions+multivariate homoscedasticity
# CONVERTED FROM MATLAB WITH THE FOLLOWING COMMENTS:
# Multivariate Statistical Testing for the Homogeneity of Covariance Matrices by the Box's M. 
# 
# Output: 
# MBox - the Box's M statistic. 
# Chi-sqr. or F - the approximation statistic test. 
# df's - degrees' of freedom of the approximation statistic test. 
# P - observed significance level. 
# 
# If the groups sample-size is at least 20 (sufficiently large), 
# Box's M test takes a Chi-square approximation; otherwise it takes 
# an F approximation. 
# 
# Example: For a two groups (g = 2) with three independent variables 
# (p = 3), we are interested in testing the homogeneity of covariances 
# matrices with a significance level = 0.05. The two groups have the 
# same sample-size n1 = n2 = 5. 
# Group 
# --------------------------------------- 
# 1 2 
# --------------------------------------- 
# x1 x2 x3 x1 x2 x3 
# --------------------------------------- 
# 23 45 15 277 230 63 
# 40 85 18 153 80 29 
# 215 307 60 306 440 105 
# 110 110 50 252 350 175 
# 65 105 24 143 205 42 
# --------------------------------------- 
# 
# Total data matrix must be: 
# X=[1 23 45 15;1 40 85 18;1 215 307 60;1 110 110 50;1 65 105 24; 
# 2 277 230 63;2 153 80 29;2 306 440 105;2 252 350 175;2 143 205 42]; 
# 
# Calling on Matlab the function: 
# MBoxtest(X,0.05) 
# 
# Answer is: 
# 
# ------------------------------------------------------------ 
# MBox F df1 df2 P 
# ------------------------------------------------------------ 
# 27.1622 2.6293 6 463 0.0162 
# ------------------------------------------------------------ 
# Covariance matrices are significantly different. 
# 
# Created by A. Trujillo-Ortiz and R. Hernandez-Walls 
# Facultad de Ciencias Marinas 
# Universidad Autonoma de Baja California 
# Apdo. Postal 453 
# Ensenada, Baja California 
# Mexico. 
# atrujo_at_uabc.mx 
# And the special collaboration of the post-graduate students of the 
# Multivariate Statistics Course: Karel Castro-Morales, 
# Alejandro Espinoza-Tenorio, Andrea Guia-Ramirez, Raquel Muniz-Salazar, 
# Jose Luis Sanchez-Osorio and Roberto Carmona-Pina. 
# November 2002. 
# 
# To cite this file, this would be an appropriate format: 
# Trujillo-Ortiz, A., R. Hernandez-Walls, K. Castro-Morales, 
# A. Espinoza-Tenorio, A. Guia-Ramirez and R. Carmona-Pina. (2002). 
# MBoxtest: Multivariate Statistical Testing for the Homogeneity of 
# Covariance Matrices by the Box's M. A MATLAB file. [WWW document]. 
# URL 
# http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=2733&objectType=FILE 
# 
# References: 
# 
# Stevens, J. (1992), Applied Multivariate Statistics for Social Sciences. 
# 2nd. ed., New-Jersey:Lawrance Erlbaum Associates Publishers. pp. 260-269. 
#----------------------------------------------------#

BoxsMTest <- function(X, vars, factorVar, alpha=0.05)
  UseMethod("BoxsMTest")
  
BoxsMTest.default <- function(X, vars, factorVar, alpha=0.05) { 
   
  x <- na.omit(X)
  if (alpha <= 0 || alpha >= 1) 
    stop('significance level must be between 0 and 1') 
  X[,factorVar] = factor(X[,factorVar])    
  g = nlevels(X[,factorVar]) ##g = nlevels(cl) ## Number of groups. 
  n = table(X[,factorVar]) ##n = table(cl) ## Vector of groups-size. 
  N = nrow(X) 
  p = length(vars) ##ncol(X) # group
  bandera = 2 
  if (any(n >= 20)) bandera = 1 
  ## Partition of the group covariance matrices. 
  ## covList <- tapply(X, rep(cl, ncol(X)), function(x, nc) cov(matrix(x, nc=nc)), ncol(X)) 

  covMat=array(dim=c(p,p,g))
  for (i in 1:g) { 
    covMat[,,i] = cov(X[X[,factorVar]==levels(X[,factorVar])[i],vars])
  }
    
  deno = sum(n) - g 
  suma = array(0, dim=dim(covMat[,,1])) ###suma = array(0, dim=dim(covList[[1]])) 
  for (k in 1:g) 
    suma = suma + (n[k] - 1) * covMat[,,k]###covList[[k]] 
  Sp = suma / deno ## Pooled covariance matrix. 
  Falta=0 
  for (k in 1:g) 
    Falta = Falta + ((n[k] - 1) * log(det(covMat[,,k]))) ###Falta = Falta + ((n[k] - 1) * log(det(covList[[k]]))) 

  MB = (sum(n) - g) * log(det(Sp)) - Falta ## Box's M statistic. 
  suma1 = sum(1 / (n[1:g] - 1)) 
  suma2 = sum(1 / ((n[1:g] - 1)^2)) 
  C = (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) * (g - 1))) * 
    (suma1 - (1 / deno)) ## Computing of correction factor. 
  if (bandera == 1) { 
    X2 = MB * (1 - C) ## Chi-square approximation. 
    v = as.integer((p * (p + 1) * (g - 1)) / 2) ## Degrees of freedom. 
    ## Significance value associated to the observed Chi-square statistic. 
    P = pchisq(X2, v, lower.tail=TRUE) 
#     cat('------------------------------------------------\n'); 
#     cat(' MBox Chi-sqr. df P\n') 
#     cat('------------------------------------------------\n') 
#     cat(sprintf("%10.4f%11.4f%12.i%13.4f\n", MB, X2, v, P)) 
#     cat('------------------------------------------------\n') 
#     if (P >= alpha) { 
#       cat('Covariance matrices are not significantly different.\n') 
#     } else { 
#       cat('Covariance matrices are significantly different.\n') 
#     } 
    return(list(MBox=MB, ChiSq=X2, df=v, pValue=P)) 
  } else { 
  ## To obtain the F approximation we first define Co, which combined to 
  ## the before C value are used to estimate the denominator degrees of 
  ## freedom (v2); resulting two possible cases. 
    Co = (((p-1) * (p+2)) / (6 * (g-1))) * (suma2 - (1 / (deno^2))) 
    if (Co - (C^2) >= 0) { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator DF. 
      v21 = as.integer(trunc((v1 + 2) / (Co - (C^2)))) ## Denominator DF. 
      F1 = MB * ((1 - C - (v1 / v21)) / v1) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P1 = pf(F1, v1, v21, lower.tail=FALSE) 
#       cat('\n------------------------------------------------------------\n') 
#       cat(' MBox F df1 df2 P\n') 
#       cat('------------------------------------------------------------\n') 
#       cat(sprintf("%10.4f%11.4f%11.i%14.i%13.4f\n", MB, F1, v1, v21, P1)) 
#       cat('------------------------------------------------------------\n') 
#       if (P1 >= alpha) { 
#         cat('Covariance matrices are not significantly different.\n') 
#       } else { 
#         cat('Covariance matrices are significantly different.\n') 
#       } 
      return(list(MBox=MB, F=F1, df1=v1, df2=v21, pValue=P1)) 
    } else { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator df. 
      v22 = as.integer(trunc((v1 + 2) / ((C^2) - Co))) ## Denominator df. 
      b = v22 / (1 - C - (2 / v22)) 
      F2 = (v22 * MB) / (v1 * (b - MB)) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P2 = pf(F2, v1, v22, lower.tail=FALSE) 
#       cat('\n------------------------------------------------------------\n') 
#       cat(' MBox F df1 df2 P\n') 
#       cat('------------------------------------------------------------\n') 
#       cat(sprintf('%10.4f%11.4f%11.i%14.i%13.4f\n', MB, F2, v1, v22, P2)) 
#       cat('------------------------------------------------------------\n') 
#       if (P2 >= alpha) { 
#         cat('Covariance matrices are not significantly different.\n') 
#       } else { 
#         cat('Covariance matrices are significantly different.\n') 
#       } 
      return(list(MBox=MB, F=F2, df1=v1, df2=v22, pValue=P2)) 
    } 
  } 
} 
