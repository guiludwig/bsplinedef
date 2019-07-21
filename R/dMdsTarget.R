dMdsTarget <- function(theta,
                       DF1 = df1, DF2 = df2,
                       M = model0, X = x, Y = hatf0, w = W){
  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1*DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1*DF2) 
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  dssq <- -2*crossprod(w, Y) + 2*crossprod(w, cbind(f1, f2))
  return(as.numeric(dssq))
}