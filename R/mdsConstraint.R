mdsConstraint <- function(theta,
                          DF1 = df1, DF2 = df2,
                          M = model0, X = x, Y = hatf0, w = W){
  # to be called as eval_g_ineq
  K2 <- length(theta)/2
  K <- sqrt(K2)
  theta1 <- matrix(theta[1:K2], ncol = K)
  theta2 <- matrix(theta[K2 + 1:K2], ncol = K)
  restrictions <- numeric(4*(K-1)*(K-1))
  jacobian <- matrix(0, nrow = 4*(K-1)*(K-1), ncol = length(theta))
  h <- 1
  for(j in 2:K){
    for(i in 2:K){
      # Restrictions
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
      # Jacobian theta 1
      jacobian[h, i - 2 + ((j-1)-2)*(K-1)] <- 0
      jacobian[h, i - 1 + ((j-1)-2)*(K-1)] <- theta2[i, j] - theta2[i-1, j]
      jacobian[h, i - 2 + ((j-1)-1)*(K-1)] <- theta2[i, j-1] - theta2[i, j]
      jacobian[h, i - 1 + ((j-1)-1)*(K-1)] <- theta2[i-1, j] - theta2[i, j-1]
      jacobian[h + 1, i - 2 + ((j-1)-2)*(K-1)] <- theta2[i-1, j] - theta2[i, j]
      jacobian[h + 1, i - 1 + ((j-1)-2)*(K-1)] <- 2*(theta2[i, j] - theta2[i-1, j])
      jacobian[h + 1, i - 2 + ((j-1)-1)*(K-1)] <- 2*theta2[i, j-1] - theta2[i, j] - theta2[i-1, j-1]
      jacobian[h + 1, i - 1 + ((j-1)-1)*(K-1)] <- -2*theta2[i, j-1] + theta2[i-1, j] + theta2[i-1, j-1] 
      jacobian[h + 2, i - 2 + ((j-1)-2)*(K-1)] <- theta2[i, j] - theta2[i, j-1]
      jacobian[h + 2, i - 1 + ((j-1)-2)*(K-1)] <- -2*theta2[i-1, j] + theta2[i, j] + theta2[i-1, j-1]
      jacobian[h + 2, i - 2 + ((j-1)-1)*(K-1)] <- -2*theta2[i, j] + 2*theta2[i, j-1]
      jacobian[h + 2, i - 1 + ((j-1)-1)*(K-1)] <- 2*theta2[i-1, j] - theta2[i, j-1] - theta2[i-1, j-1]
      jacobian[h + 3, i - 2 + ((j-1)-2)*(K-1)] <- theta2[i-1, j] - theta2[i, j-1]
      jacobian[h + 3, i - 1 + ((j-1)-2)*(K-1)] <- 2*(theta2[i, j] - theta2[i-1, j]) - theta2[i-1, j] + theta2[i-1, j-1]
      jacobian[h + 3, i - 2 + ((j-1)-1)*(K-1)] <- 3*theta2[i, j-1] - 2*theta2[i, j] - theta2[i-1, j-1] 
      jacobian[h + 3, i - 1 + ((j-1)-1)*(K-1)] <- 2*theta2[i-1, j] - 2*theta2[i, j-1]
      # Jacobian theta 2
      jacobian[h, (K-1)^2 + i - 2 + ((j-1)-2)*(K-1)] <- 0
      jacobian[h, (K-1)^2 + i - 1 + ((j-1)-2)*(K-1)] <- theta1[i-1, j] - theta1[i, j]
      jacobian[h, (K-1)^2 + i - 2 + ((j-1)-1)*(K-1)] <- theta1[i, j] - theta1[i, j-1]
      jacobian[h, (K-1)^2 + i - 1 + ((j-1)-1)*(K-1)] <- theta1[i, j-1] - theta1[i-1, j]
      jacobian[h + 1, (K-1)^2 + i - 2 + ((j-1)-2)*(K-1)] <- theta1[i,j] - theta1[i-1,j]
      jacobian[h + 1, (K-1)^2 + i - 1 + ((j-1)-2)*(K-1)] <- 2*(theta1[i-1,j] - theta1[i,j])
      jacobian[h + 1, (K-1)^2 + i - 2 + ((j-1)-1)*(K-1)] <- -2*theta1[i,j-1] + theta1[i-1,j-1] + theta1[i-1,j]
      jacobian[h + 1, (K-1)^2 + i - 1 + ((j-1)-1)*(K-1)] <- 2*theta1[i,j-1] - theta1[i-1,j] - theta1[i-1,j-1]
      jacobian[h + 2, (K-1)^2 + i - 2 + ((j-1)-2)*(K-1)] <- theta1[i, j-1] - theta1[i, j]
      jacobian[h + 2, (K-1)^2 + i - 1 + ((j-1)-2)*(K-1)] <- 2*theta1[i-1, j] - theta1[i, j] - theta1[i-1, j-1]
      jacobian[h + 2, (K-1)^2 + i - 2 + ((j-1)-1)*(K-1)] <- -2*theta1[i, j-1] + 2*theta1[i, j]
      jacobian[h + 2, (K-1)^2 + i - 1 + ((j-1)-1)*(K-1)] <- -2*theta1[i-1, j] + theta1[i, j-1] + theta1[i-1, j-1]
      jacobian[h + 3, (K-1)^2 + i - 2 + ((j-1)-2)*(K-1)] <- -theta1[i-1, j] + theta1[i, j-1]
      jacobian[h + 3, (K-1)^2 + i - 1 + ((j-1)-2)*(K-1)] <- -2*(theta1[i, j] - theta1[i-1, j]) + theta1[i-1, j] - theta1[i-1, j-1]
      jacobian[h + 3, (K-1)^2 + i - 2 + ((j-1)-1)*(K-1)] <- -3*theta1[i, j-1] + 2*theta1[i, j] + theta1[i-1, j-1] 
      jacobian[h + 3, (K-1)^2 + i - 1 + ((j-1)-1)*(K-1)] <- -2*theta1[i-1, j] + 2*theta1[i, j-1]
      # Increase counter
      h <- h + 4
    }
  }
  return(list(constraints = -1*restrictions,
              jacobian = -1*jacobian))
  # return(-1*restrictions)
}