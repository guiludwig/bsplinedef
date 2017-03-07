likelihoodTarget <- function(theta,
                             DF1 = df1, DF2 = df2, B = basis,
                             M = model0, X = x, Y = y){

  n <- nrow(X) #!# length(Y)
  W <- matrix(0, n, DF1*DF2)
  # Need to turn around!
  for(i in 1:n) W[i,] <- kronecker(predict(B$B2, X[, 2])[i, ],
                                   predict(B$B1, X[, 1])[i, ])
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  f1 <- as.numeric(W%*%theta1)
  f2 <- as.numeric(W%*%theta2)
  lik <- -1*RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood
  return(lik)
}
