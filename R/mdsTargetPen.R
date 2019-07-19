mdsTargetPen <- function(theta,
                         DF1 = df1, DF2 = df2,
                         M = model0, X = x, Y = hatf0, w = W){
  n <- nrow(X)
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  likP <- sum((f1 - Y[,1])^2) + sum((f2 - Y[,2])^2) 
  return(likP)
}