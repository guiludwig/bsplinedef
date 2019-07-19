mdsConstraint <- function(theta,
                          DF1 = df1, DF2 = df2,
                          M = model0, X = x, Y = hatf0, w = W){
  # to be called as eval_g_ineq
  K2 <- length(theta)/2
  K <- sqrt(K2)
  theta1 <- matrix(theta[1:K2], ncol = K)
  theta2 <- matrix(theta[K2 + 1:K2], ncol = K)
  restrictions <- numeric(4*(K-1)*(K-1))
  h <- 1
  for(i in 2:K){
    for(j in 2:K){
      restrictions[h] <- theta1[i,j-1]*(theta2[i, j]-theta2[i-1, j]) - 
        theta1[i,j]*(theta2[i,j-1]-theta2[i-1, j-1]) +
        theta1[i,j]*(theta2[i-1,j]-theta2[i-1, j-1]) -
        theta1[i-1,j]*(theta2[i,j]-theta2[i, j-1])
      restrictions[h+1] <- restrictions[h] +
        theta1[i,j-1]*(theta2[i, j]-theta2[i-1, j]) - 
        theta1[i,j]*(theta2[i,j-1]-theta2[i-1, j-1]) +
        theta1[i-1,j]*(theta2[i, j-1]-theta2[i-1, j-1]) - 
        theta1[i-1,j-1]*(theta2[i,j]-theta2[i-1, j])
      restrictions[h+2] <- restrictions[h] +
        theta1[i,j]*(theta2[i-1,j]-theta2[i-1, j-1]) -
        theta1[i-1,j]*(theta2[i,j]-theta2[i, j-1]) + 
        theta1[i-1,j-1]*(theta2[i,j]-theta2[i, j-1]) -
        theta1[i,j-1]*(theta2[i-1,j]-theta2[i-1, j-1])
      restrictions[h+3] <- 2*restrictions[h] + 
        theta1[i-1,j]*(theta2[i, j-1]-theta2[i-1, j-1]) - 
        theta1[i-1,j-1]*(theta2[i,j]-theta2[i-1, j]) +
        theta1[i-1,j-1]*(theta2[i,j]-theta2[i, j-1]) -
        theta1[i,j-1]*(theta2[i-1,j]-theta2[i-1, j-1])
      h <- h + 4
    }
  }
  # return(list(constraints = -1*restrictions,
  #             jacobian = ??))
  return(-1*restrictions)
}