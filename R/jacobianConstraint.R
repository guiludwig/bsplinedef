jacobianConstraint <- function(theta,
                               DF1 = df1, DF2 = df2, B = basis, # After here not used, match?
                               b1 = B1, b2 = B2, db1 = dB1, db2 = dB2,
                               M = model0, X = x, Y = y, ...){

  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  PEN <- numeric(nrow(X))

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

  for(i in 1:n){
    df1d1 <- db1[i,]%*%theta1%*%b2[i,]
    df2d2 <- b1[i,]%*%theta2%*%db2[i,]
    df1d2 <- b1[i,]%*%theta1%*%db2[i,]
    df2d1 <- db1[i,]%*%theta2%*%b2[i,]
    PEN[i] <- as.numeric(df1d1*df2d2 - df1d2*df2d1)
  }

  return(PEN)

}
