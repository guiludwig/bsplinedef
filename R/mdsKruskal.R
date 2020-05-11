mdsKruskal <- function(x0, sv = SV, M = model1){
  x0 <- matrix(x0, ncol = 2) # will be passed as numeric
  a <- summary(M)
  sv <- sv/2 # Semivariogram
  # Missing ROOT!
  if(any(grep("nugget", colnames(a$param)))){
    sv[lower.tri(sv)] <- sv[lower.tri(sv)] - a$param["value", grep("nugget", colnames(a$param))]
    sv[upper.tri(sv)] <- sv[upper.tri(sv)] - a$param["value", grep("nugget", colnames(a$param))]
    sv[sv < 0] <- 0
  }
  a.nuggetless <- a$param[,-grep("nugget", colnames(a$param))]
  if(length(dim(a.nuggetless)) <= 1){
    ret <- MASS::isoMDS(sv, trace = FALSE)$points
  } else {
    # Assuming largest spatial effect is most relevant
    temp <- a.nuggetless["value", grep("\\.var$", colnames(a.nuggetless))]
    temp2 <- a.nuggetless["value",grep("\\.s$", colnames(a.nuggetless))]
    sigma2 <- max(temp)
    sv[sv > 0.99*sigma2] <- 0.99*sigma2 # Small proportion, 0.008 in simulation
    phi <- temp2[order(temp)][1]
    sv <- 1 - sv/sigma2
    dij <- -phi*log(sv)
    ret <- cmdscale(dij)
  }
  # Cdelta <- RandomFields::RFcovmatrix(M, x0)
  # GammaDelta <- kronecker(matrix(1, nrow = nrow(Cdelta)),
  #                         matrix(diag(Cdelta), ncol = ncol(Cdelta))) +
  #   kronecker(matrix(1, ncol = ncol(Cdelta)),
  #             matrix(diag(Cdelta), nrow = nrow(Cdelta))) -2*Cdelta
  # if(any(dim(GammaDelta) != dim(sv))) stop("Matrix dimension mismatch.")
  # Sampson and Guttorp's 1/|x_i - x_j| weight
  # WEIGHT <- as.matrix(dist(x0))
  # diag(WEIGHT) <- rep(1, nrow(WEIGHT))
  # WEIGHT <- 1/WEIGHT
  # diag(WEIGHT) <- rep(0, nrow(WEIGHT))
  # ssq <- sum(WEIGHT*(GammaDelta - sv)^2)/sum(WEIGHT*(GammaDelta)^2)
  # ssq <- sum((GammaDelta - sv)^2)/sum((GammaDelta)^2)
  return(ret)
}