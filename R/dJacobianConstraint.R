dJacobianConstraint <- function(theta,
                                DF1 = df1, DF2 = df2, B = basis, # After here not used, match?
                                M = model0, X = x, Y = y){

  # Since jacobianConstraint is f : R^K -> R^n,
  # this needs to return a n x K matrix.
  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  dPEN <- matrix(nrow = nrow(X), ncol = 2*DF1*DF2)

  dB1x1 <- splineDesign(knots = c(rep(attr(B$B1, "Boundary.knots")[1],
                                      attr(B$B1, "degree")+1),
                                  attr(B$B1, "knots"),
                                  rep(attr(B$B1, "Boundary.knots")[2],
                                      attr(B$B1, "degree")+1)),
                        x = X[, 1, drop = FALSE], outer.ok = TRUE,
                        derivs = 1)
  dB2x2 <- splineDesign(knots = c(rep(attr(B$B2, "Boundary.knots")[1],
                                      attr(B$B2, "degree")+1),
                                  attr(B$B2, "knots"),
                                  rep(attr(B$B2, "Boundary.knots")[2],
                                      attr(B$B2, "degree")+1)),
                        x = X[, 2, drop = FALSE], outer.ok = TRUE,
                        derivs = 1)
  B1x1 <- predict(B$B1, X[ , 1, drop = FALSE])
  B2x2 <- predict(B$B2, X[ , 2, drop = FALSE])

  for(i in 1:n){
    for(k1 in 1:DF1){
      for(k2 in 1:DF2){
        dPEN[i, k1 + (k2-1)*DF1] <- B1x1[i,k1]*B2x2[i,k2]*(dB1x1[i,k1]*crossprod(dB2x2[i,], theta2[k1,]) - dB2x2[i,k2]*crossprod(dB1x1[i,], theta2[,k2]))
        dPEN[i, DF1*DF2 + k1 + (k2-1)*DF1] <- -1*B1x1[i,k1]*B2x2[i,k2]*(dB1x1[i,k1]*crossprod(dB2x2[i,], theta1[k1,]) - dB2x2[i,k2]*crossprod(dB1x1[i,], theta1[,k2]))
      }
    }
  }

  return(dPEN)

}
