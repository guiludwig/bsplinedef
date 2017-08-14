dJacobianConstraint <- function(theta,
                                DF1 = df1, DF2 = df2, B = basis, # After here not used, match?
                                b1 = B1, b2 = B2, db1 = dB1, db2 = dB2,
                                M = model0, X = x, Y = y, ...){

  # Since jacobianConstraint is f : R^K -> R^n,
  # this needs to return a n x K matrix.
  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  dPEN <- matrix(nrow = nrow(X), ncol = 2*DF1*DF2)

  # dB1 <- splineDesign(knots = c(rep(attr(B$B1, "Boundary.knots")[1],
  #                                   attr(B$B1, "degree")+1),
  #                               attr(B$B1, "knots"),
  #                               rep(attr(B$B1, "Boundary.knots")[2],
  #                                   attr(B$B1, "degree")+1)),
  #                     x = X[, 1, drop = FALSE], outer.ok = TRUE,
  #                     derivs = 1)
  # dB2 <- splineDesign(knots = c(rep(attr(B$B2, "Boundary.knots")[1],
  #                                   attr(B$B2, "degree")+1),
  #                               attr(B$B2, "knots"),
  #                               rep(attr(B$B2, "Boundary.knots")[2],
  #                                   attr(B$B2, "degree")+1)),
  #                     x = X[, 2, drop = FALSE], outer.ok = TRUE,
  #                     derivs = 1)
  # B1 <- predict(B$B1, X[ , 1, drop = FALSE])
  # B2 <- predict(B$B2, X[ , 2, drop = FALSE])

  for(ell in 1:n){
    for(i in 1:DF1){
      for(j in 1:DF2){
        dPEN[ell, i + (j-1)*DF1] <- b1[ell,i]*b2[ell,j]*(db1[ell,i]*crossprod(db2[ell,], theta2[i,]) -
                                                           db2[ell,j]*crossprod(db1[ell,], theta2[,j]))
        dPEN[ell, DF1*DF2 + i + (j-1)*DF1] <- -1*b1[ell,i]*b2[ell,j]*(db1[ell,i]*crossprod(db2[ell,], theta1[i,]) -
                                                                        db2[ell,j]*crossprod(db1[ell,], theta1[,j]))
      }
    }
  }

  return(dPEN)

}
