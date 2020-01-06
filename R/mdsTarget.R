mdsTarget <- function(theta,
                      DF1 = df1, DF2 = df2,
                      M = model0, X = x, Y = hatf0, w = W){
  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1*DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1*DF2) 
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  resid <- Y - cbind(f1,f2)
  # Frobenius norm
  ssq <- sum(diag(crossprod(resid)))
  return(ssq)
}