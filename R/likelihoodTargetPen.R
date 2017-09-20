likelihoodTargetPen <- function(theta,
                                DF1 = df1, DF2 = df2,
                                M = model0, X = x, Y = y, w = W, 
                                LAMBDA = lambda,  
                                z1 = zeta1, z2 = zeta2){
  n <- nrow(X)
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  likP <- -1*RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood + 
    LAMBDA * injectSuff(matrix(theta1, nrow = DF1, ncol = DF2), 
                        matrix(theta2, nrow = DF1, ncol = DF2), 
                        z1, z2)
  return(likP)
}