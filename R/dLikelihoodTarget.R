# Derivatives WITH RESPECT TO POSITION!!
dLikelihoodTarget <- function(theta,
                             DF1 = df1, DF2 = df2,
                             M = model0, X = x, Y = y, w = W, ...){

  # Note d/d(theta1) log-L(f1,f2) = d/df log-L(f1,f2) df/dtheta
  # RFcov?
  n <- nrow(X) #!# length(Y)
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  lik <- -1*RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood
  return(lik)
}
