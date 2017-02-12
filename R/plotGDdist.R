#' Scatterplot of sample spatial covariance entries and pairwise distances in G- and D-domain
#'
#' This function produces a plot of the upper diagonal sample covariance entries and the
#' corresponding pairwise distances in the original and deformed coordinates, where the
#' spatial deformation is obtained with a tensor product of B-spline basis functions,
#' implemented in the \code{\link{bdef}} function.
#'
#' @export
plotGDdist <- function(model, plot = TRUE, lines = TRUE) {

  x <- model$x
  y <- model$y
  spatial.locations <- unique(x)
  p <- nrow(spatial.locations)
  fine <- (length(p) > 0)&&(p>2)
  if(!fine) stop("Error evaluating sample covariance matrix. Make sure that there are more than one observation at each sampled site.")
  sample.sigma <- matrix(0, p, p)
  for(i in 1:p){
    for(j in i:p){
      index.i <- which(x[,1] == spatial.locations[i,1] & x[,2] == spatial.locations[i,2])
      index.j <- which(x[,1] == spatial.locations[j,1] & x[,2] == spatial.locations[j,2])
      if(length(index.i) == 1 | length(index.j) == 1) stop("Sites must have more than one observation.")
      tryCatch(sample.sigma[i,j] <- cov(y[index.i], y[index.j]),
               error = function(e) simpleError(message = "sites must have equal number of observations."))
    }
  }
  sample.sigma <- sample.sigma + t(sample.sigma)

  cov.model <- model$model
  Gcoords <- dist(model$x)
  Gcoords <- as.matrix(Gcoords)
  Dcoords <- dist(model$def.x)
  Dcoords <- as.matrix(Dcoords)
  # G <- RFcovmatrix(model, distances = Gcoords[lower.tri(Gcoords)], dim = nrow(Gcoords))
  # D <- RFcovmatrix(model, distances = Dcoords[lower.tri(Dcoords)], dim = nrow(Dcoords))
  layout(matrix(1:2, ncol = 2))
  plot(as.numeric(Gcoords), sample.sigma[lower.tri(sample.sigma)],
       main = "G-domain")
  plot(as.numeric(Dcoords), sample.sigma[lower.tri(sample.sigma)],
       main = "D-domain")
  # TO IMPLEMENT: LINES

}
