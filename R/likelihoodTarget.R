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
  if(is.na(L) & length(theta) == 2*DF1*DF2 + 1){
    L <- theta[2*DF1*DF2+1]
  }
  f1 <- as.numeric(W%*%theta1)
  f2 <- as.numeric(W%*%theta2)
  # lik <- lik + lambda*PEN
  dB1 <- splineDesign(knots = c(rep(attr(B$B1, "Boundary.knots")[1],
                                    attr(B$B1, "degree")+1),
                                attr(B$B1, "knots"),
                                rep(attr(B$B1, "Boundary.knots")[2],
                                    attr(B$B1, "degree")+1)),
                      x = X[, 1], outer.ok = TRUE,
                      derivs = 1)
  dB2 <- splineDesign(knots = c(rep(attr(B$B2, "Boundary.knots")[1],
                                    attr(B$B2, "degree")+1),
                                attr(B$B2, "knots"),
                                rep(attr(B$B2, "Boundary.knots")[2],
                                    attr(B$B2, "degree")+1)),
                      x = X[, 2], outer.ok = TRUE,
                      derivs = 1)
  dWdx1 <- matrix(0, n, DF1*DF2)
  for(i in 1:n) dWdx1[i,] <- kronecker(dB1[i, ],
                                       predict(B$B2, X[, 2])[i, ])
  dWdx2 <- matrix(0, n, DF1*DF2)
  for(i in 1:n) dWdx2[i,] <- kronecker(predict(B$B1, X[, 1])[i, ],
                                       dB2[i, ])
  df1dx1 <- as.numeric(dWdx1%*%theta1)
  df1dx2 <- as.numeric(dWdx2%*%theta1)
  df2dx1 <- as.numeric(dWdx1%*%theta2)
  df2dx2 <- as.numeric(dWdx2%*%theta2)
  df1dx1 <- scale(df1dx1, center = FALSE, scale = TRUE)
  df1dx2 <- scale(df1dx2, center = FALSE, scale = TRUE)
  df2dx1 <- scale(df2dx1, center = FALSE, scale = TRUE)
  df2dx2 <- scale(df2dx2, center = FALSE, scale = TRUE)
  # TODO: debug: is this right?
  PEN <- abs(sum(df1dx1*df2dx2) - sum(df1dx2*df1dx2))
  lik <- -1*RFlikelihood(M, x = f1, y = f2, data = Y)$loglikelihood + exp(L)*PEN
  return(lik)
}
