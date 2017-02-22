jacobianConstraint <- function(theta,
         DF1 = df1, DF2 = df2,
         B = basis, X = x, ...){

  n <- nrow(X)
  theta1 <- matrix(theta[1:(DF1*DF2)], nrow = DF1, ncol = DF2)
  theta2 <- matrix(theta[1:(DF1*DF2) + (DF1*DF2)], nrow = DF1, ncol = DF2)
  PEN <- numeric(nrow(X))

  dB1 <- splineDesign(knots = c(rep(attr(B$B1, "Boundary.knots")[1],
                                    attr(B$B1, "degree")+1),
                                attr(B$B1, "knots"),
                                rep(attr(B$B1, "Boundary.knots")[2],
                                    attr(B$B1, "degree")+1)),
                      x = X[, 1, drop = FALSE], outer.ok = TRUE,
                      derivs = 1)
  dB2 <- splineDesign(knots = c(rep(attr(B$B2, "Boundary.knots")[1],
                                    attr(B$B2, "degree")+1),
                                attr(B$B2, "knots"),
                                rep(attr(B$B2, "Boundary.knots")[2],
                                    attr(B$B2, "degree")+1)),
                      x = X[, 2, drop = FALSE], outer.ok = TRUE,
                      derivs = 1)
  B1 <- predict(B$B1, X[ , 1, drop = FALSE])
  B2 <- predict(B$B2, X[ , 2, drop = FALSE])

  for(i in 1:n){
    t1b1 <- theta1%*%B1[i,]
    t1b2 <- theta1%*%B2[i,]
    t2b1 <- theta2%*%B1[i,]
    t2b2 <- theta2%*%B2[i,]
    PEN[i] <- as.numeric(B1[i,]%*%(tcrossprod(t2b2,t1b1) - tcrossprod(t1b2,t2b1))%*%B2[i,])
  }

  return(PEN)

}
