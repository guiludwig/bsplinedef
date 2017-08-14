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
    # t1b1 <- theta1%*%dB1[i,]
    # t1b2 <- theta1%*%dB2[i,]
    # t2b1 <- theta2%*%dB1[i,]
    # t2b2 <- theta2%*%dB2[i,]
    df1d1 <- db1[i,]%*%theta1%*%b2[i,]
    df2d2 <- b1[i,]%*%theta2%*%db2[i,]
    df1d2 <- b1[i,]%*%theta1%*%db2[i,]
    df2d1 <- db1[i,]%*%theta2%*%b2[i,]
    # Test: both eigenvalues have to be positive
    TT <- as.numeric(df1d1 + df2d2)
    DD <- as.numeric(df1d1*df2d2 - df1d2*df2d1)
    L1 <- TT/2 + sqrt(as.complex((TT/2)^2-DD))
    L2 <- TT/2 - sqrt(as.complex((TT/2)^2+DD))
    #!# PEN[i] <- as.numeric(df1d1*df2d2 - df1d2*df2d1) # Same as
    if(Im(L1) != 0 | Im(L2) != 0) {
      PEN[i] <- -abs(Re(L1)*Re(L2)) - abs(Im(L1)*Im(L2)) # DD
    } else {
      PEN[i] <- Re(L1)*Re(L2) # DD
    }
    # Profile later
    # PEN[i] <- as.numeric(B1[i,]%*%(theta2%*%outer(dB2[i,],dB1[i,])%*%theta1 -
    #                                theta1%*%outer(dB2[i,],dB1[i,])%*%theta2)%*%B2[i,])
    # PEN[i] <- as.numeric(B1[i,]%*%(tcrossprod(t2b2,t1b1) - tcrossprod(t1b2,t2b1))%*%B2[i,])
  }

  return(PEN)

}
