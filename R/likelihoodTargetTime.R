likelihoodTargetTime <- function(theta,
                             DF1 = df1, DF2 = df2, B = basis,
                             M = model0, X = x, Y = y, TIM = tim){

  # For now: treat reps as independent samples
  n <- nrow(X) #!# length(Y)
  m <- length(TIM)
  W <- matrix(0, n, DF1*DF2)
  #!# for(i in 1:n) W[i,] <- kronecker(predict(B$B1, rep(X[, 1], m))[i, ],
  #!#                                  predict(B$B2, rep(X[, 2], m))[i, ])
  for(i in 1:n) W[i,] <- kronecker(predict(B$B1, X[, 1])[i, ],
                                   predict(B$B2, X[, 2])[i, ])
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  # if(is.na(L) & length(theta) == 2*DF1*DF2 + 1){
  #   L <- theta[2*DF1*DF2+1]
  # }
  f1 <- as.numeric(W%*%theta1)
  f2 <- as.numeric(W%*%theta2)
  lik <- -1*RFlikelihood(M, x = f1, y = f2, data = matrix(Y, ncol = m))$loglikelihood # + exp(L)*PEN
  return(lik)
}
