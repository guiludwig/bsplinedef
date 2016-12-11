likelihoodTarget <- function(theta, DF1 = df1, DF2 = df2, B = basis,
                             M = model0, X = x, Y = y, L = lambda){
  # coords <- dist(X)
  # coords <- as.matrix(coords)
  # # Evals cov:
  # SIGMA <- RFcovmatrix(model0, distances = coords[lower.tri(coords)], dim = nrow(coords))
  # SIGMA <- (SIGMA + t(SIGMA))/2 # Ensures symmetry after rounding error? Small dim is fine
  n <- length(Y)
  W <- matrix(0, n, DF1*DF2)
  for(i in 1:n) W[i,] <- kronecker(predict(B$B1, X[, 1])[i, ],
                                   predict(B$B2, X[, 2])[i, ])
  theta1 <- theta[1:(DF1*DF2)]
  theta2 <- theta[1:(DF1*DF2) + (DF1*DF2)]
  if(is.na(L)){
    L <- theta[2*DF1*DF2+1]
  }
  f1 <- as.numeric(W%*%theta1)
  f2 <- as.numeric(W%*%theta2)
  # lik <- lik + lambda*PEN
  lik <- RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood + L*1
  #
  # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE,
  #              sparse = FALSE)
  return(lik)
}
