predict.bdef <- function(object, newdata) {
  # TODO: check if newdata in object$window and/or matrix
  b1 <- predict(object$basis$B1, newdata[,1])
  b2 <- predict(object$basis$B2, newdata[,2])
  # Retrive coefficients
  theta1 <- matrix(object$theta1, nrow = object$df1, ncol = object$df2)
  theta2 <- matrix(object$theta2, nrow = object$df1, ncol = object$df2)
  # Deform new coordinates
  nf1 <- nf2 <- numeric(nrow(newdata))
    for(i in 1:nrow(newdata)){
    nf1[i] <- b1[i,]%*%theta1%*%b2[i,]
    nf2[i] <- b1[i,]%*%theta2%*%b2[i,]
  }
  newdata.def <- round(cbind(nf1, nf2),6)
  # Predict here
  # RFinterpolate(object$model, x=nf1, y=nf2) # UGH
  tempX <- unique(rbind(object$def.x, newdata.def)) # Because of error in repeated places
  S <- RFcovmatrix(object$model, x = tempX[,1], y = tempX[,2])
  dataRows <- which((tempX[,1] %in% object$def.x[,1]) & (tempX[,2] %in% object$def.x[,2]))
  newDataRows <- which((tempX[,1] %in% newdata.def[,1]) & (tempX[,2] %in% newdata.def[,2]))
  Sigma <- S[dataRows, dataRows]
  SigmaCross <- S[dataRows,newDataRows]
  # Average residuals in time
  R <- apply(resid(object$model)$residuals, 1, mean)
  krige <- crossprod(SigmaCross,solve(Sigma,R))
  return(list(krige = krige, newdata.def = newdata.def))
}