mdsKruskal <- function(x0, sv = SV, M = model1){
  x0 <- matrix(x0, ncol = 2) # will be passed as numeric
  Cdelta <- RandomFields::RFcovmatrix(M, x0)
  GammaDelta <- kronecker(matrix(1, nrow = nrow(Cdelta)), 
                          matrix(diag(Cdelta), ncol = ncol(Cdelta))) + 
    kronecker(matrix(1, ncol = ncol(Cdelta)), 
              matrix(diag(Cdelta), nrow = nrow(Cdelta))) -2*Cdelta
  if(any(dim(GammaDelta) != dim(sv))) stop("Matrix dimension mismatch.")
  # Sampson and Guttorp's 1/|x_i - x_j| weight
  # WEIGHT <- as.matrix(dist(x0))
  # diag(WEIGHT) <- rep(1, nrow(WEIGHT))
  # WEIGHT <- 1/WEIGHT
  # diag(WEIGHT) <- rep(0, nrow(WEIGHT))
  # ssq <- sum(WEIGHT*(GammaDelta - sv)^2)/sum(WEIGHT*(GammaDelta)^2)
  ssq <- sum((GammaDelta - sv)^2)/sum((GammaDelta)^2)
  return(sqrt(ssq))
}