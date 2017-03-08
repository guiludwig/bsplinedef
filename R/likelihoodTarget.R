likelihoodTarget <- function(theta,
                             DF1 = df1, DF2 = df2,
                             M = model0, X = x, Y = y, w = W, ...){

  n <- nrow(X) #!# length(Y)
  # W <- matrix(0, n, DF1*DF2)
  # Need to turn around!
  #for(i in 1:n) W[i,] <- kronecker(predict(B$B2, X[, 2])[i, ],
  #                                  predict(B$B1, X[, 1])[i, ])
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  # B1 <- predict(B$B1, x[ , 1, drop = FALSE])
  # B2 <- predict(B$B2, x[ , 2, drop = FALSE])
  # theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  # theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  # f2 <- f1 <- numeric(n)
  # for(i in 1:n){
    # f1[i] <- (B1[i,])%*%theta1%*%(B2[i,])
    # f2[i] <- (B1[i,])%*%theta2%*%(B2[i,])
  # }
  # f1 <- as.numeric(f1)
  # f2 <- as.numeric(f2)
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  lik <- -1*RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood
  return(lik)
}
