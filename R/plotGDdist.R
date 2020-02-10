#' Scatterplot of sample spatial covariance entries and pairwise distances in G- and D-domain
#'
#' This function produces a plot of the upper diagonal sample covariance entries and the
#' corresponding pairwise distances in the original and deformed coordinates, where the
#' spatial deformation is obtained with a tensor product of B-spline basis functions,
#' implemented in the \code{\link{bdef}} function.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 40
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x <- cbind(x1,x2)
#' F1 <- function(x1,x2) {
#'   x1 # (x1 - 2*x2)^2/sqrt(x1^2 + 4*x2^2)
#' }
#' F2 <- function(x1,x2) {
#'   sqrt(x1*x2) # (3*x1 + x2)^2/sqrt(9*x1^2 + x2^2)
#' }
#' plotGrid(list(window = list(x = range(x1), y = range(x2)),
#'               x = x, def.x = cbind(F1(x1,x2), F2(x1,x2))),
#'          F1 = F1, F2 = F2, margins = TRUE)
#' TIME <- 20
#' covModel2 <- RMexp(proj = "space", var = 1, scale = .25) *
#'   RMnugget(proj = "time", var = 1)
#' data2 <- RFsimulate(covModel2, x = F1(x1,x2), y = F2(x1,x2), T = 1:TIME)
#' y <- as.numeric(unlist(data2@data))
#' \dontrun{
#' matplot(1:TIME, t(matrix(y, ncol = TIME)), type = "b", lty = 1, pch = 1)
#' }
#' covModel2m <- RMexp(var = NA, scale = NA) # Only space?
#' testT <- bdef(x, y, tim = 1:TIME, cov.model = covModel2m, maxit = 10) # lambda estimated
#' plotGrid(testT, margins = TRUE)
#' plotGDdist(testT)
plotGDdist <- function(model, plot = TRUE, ggplot = FALSE, lines = TRUE) {

  if(is.null(model$fullDes)) stop("Please enable the fullDes option in bdef().")
  x <- model$fullDes$x
  y <- model$fullDes$y
  spatial.locations <- unique(x)
  p <- nrow(spatial.locations)
  fine <- (length(p) > 0)&&(p>2)
  if(!fine) stop("Error evaluating sample covariance matrix. Make sure that there are more than one
                 observation at each sampled site.")
  sample.sigma <- matrix(0, p, p)
  for(i in 1:p){
    for(j in i:p){
      index.i <- which(x[,1] == spatial.locations[i,1] & x[,2] == spatial.locations[i,2])
      index.j <- which(x[,1] == spatial.locations[j,1] & x[,2] == spatial.locations[j,2])
      if(length(index.i) == 1 | length(index.j) == 1) stop("Sites must have more than one
                                                           observation.")
      tryCatch(sample.sigma[i,j] <- cov(y[index.i], y[index.j]),
               error = function(e) simpleError(message = "sites must have equal number of
                                               observations."))
    }
  }
  S <- diag(sample.sigma)
  sample.sigma <- sample.sigma + t(sample.sigma)
  diag(sample.sigma) <- S

  cov.model <- model$model
  Gcoords <- dist(model$x)
  Gcoords <- as.matrix(Gcoords)
  Dcoords <- dist(model$def.x)
  Dcoords <- as.matrix(Dcoords)
  # G <- RFcovmatrix(model, distances = Gcoords[lower.tri(Gcoords)], dim = nrow(Gcoords))
  # D <- RFcovmatrix(model, distances = Dcoords[lower.tri(Dcoords)], dim = nrow(Dcoords))
  if(ggplot){
    # Need to return cov.model?
    # see ?print.RFfit
    # predict lines!
  } else {
    layout(matrix(c(1,0,2:3), ncol = 2))
    plot(as.numeric(Gcoords)[as.numeric(Gcoords)>0],
         as.numeric(sample.sigma)[as.numeric(Gcoords)>0], # sample.sigma[lower.tri(sample.sigma)],
         main = "G-domain",
         ylab = expression("Cov("*Y(s[i])*","*Y(s[j])*")"),
         xlab = expression("|"*x[i] - x[j]*"|"))
    # TO IMPLEMENT: LINES
    plot(as.numeric(Gcoords)[as.numeric(Gcoords)>0],
         as.numeric(Dcoords)[as.numeric(Dcoords)>0],
         xlab = expression("|"*x[i] - x[j]*"|"),
         ylab = expression("|"*f[i] - f[j]*"|"))
    abline(a=0,b=1)
    plot(as.numeric(Dcoords)[as.numeric(Dcoords)>0],
         as.numeric(sample.sigma)[as.numeric(Dcoords) >0], # sample.sigma[lower.tri(sample.sigma)],
         main = "D-domain",
         ylab = expression("Cov("*Y(s[i])*","*Y(s[j])*")"),
         xlab = expression("|"*f[i] - f[j]*"|"))
  }

}
