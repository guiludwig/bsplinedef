jacobianConstraint2 <- function(theta, zeta1 = .5, zeta2 = .5, 
                                DF1 = df1, DF2 = df2, ...){
  
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  
  PEN <- sum(rhoPen(apply(theta1, 1, diff), zeta1, zeta2)) + # Differences by row
    sum(rhoPen(apply(theta1, 2, diff), zeta1, zeta2)) + # Differences by column
    sum(rhoPen(apply(theta2, 1, diff), zeta1, zeta2)) + # Differences by row
    sum(rhoPen(apply(theta2, 2, diff), zeta1, zeta2)) # Differences by column

  return(PEN)

}
