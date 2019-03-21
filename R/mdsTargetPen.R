mdsTargetPen <- function(theta,
                         DF1 = df1, DF2 = df2,
                         M = model0, X = x, Y = hatf0, w = W, 
                         LAMBDA = lambda,  
                         z1 = zeta1, z2 = zeta2){
  n <- nrow(X)
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  f1 <- as.numeric(w%*%theta1)
  f2 <- as.numeric(w%*%theta2)
  likP <- sum((f1-hatf0[,1])^2) + sum((f2-hatf0[,2])^2) + 
    LAMBDA * injectSuff(matrix(theta1, nrow = DF1, ncol = DF2), 
                        matrix(theta2, nrow = DF1, ncol = DF2), 
                        z1, z2) # bottleneck
  return(likP)
}