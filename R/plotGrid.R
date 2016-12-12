#' Plot Spatial Deformations from bdef
#'
#' This function plot a spatial deformation on a grid obtained via the bdef
#' function
#'
#' @param model An object returned by the function \code{\link{bdef}}.
#'                     Alternatively, a list with named components matching
#'                     those returned by \code{\link{bdef}}.
#' @param nx Number of ticks in x-axis for grid. Defaults to 10.
#' @param ny Number of ticks in y-axis for grid. Defaults to 10.
#' @param colorO Color of original data points and grid. Defaults to Black.
#' @param colorD Color of deformed data points and grid. Defaults to Red.
#' @param plot Whether a plot should be produced or not. Defaults to TRUE.
#' @param margins Whether a layout should be produced, with f1, f2 estimates
#'                     shown along the corresponding axes. Defaults to FALSE.
#' @param F1 Allows the user to supply a function that generates the deformed
#'                     coordinate y1. Defaults to NULL, in which case the
#'                     tensor product of B-splines is used instead.
#' @param F2 Allows the user to supply a function that generates the deformed
#'                     coordinate y2. Defaults to NULL, in which case the
#'                     tensor product of B-splines is used instead.
#'
#' @export
#' @return \code{plotGrid} invisibly returns the values of the deformation on
#'                     the grid.
#'
#' @examples
#' # Example using artificially generated data
#' set.seed(1)
#' n <- 10
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x <- cbind(x1,x2)
#' \dontrun{
#'    plot(x, xlim = c(0,1), ylim = c(0,1), pch = letters[1:n])
#' }
#' K <- 8
#' theta1 <- exp(-0.5*scale(kronecker(1:K - K/2,
#'               rep(1,K) + 5*((1:K - K/2)^2)/K^2),
#'               center = FALSE, scale = TRUE))
#' theta2 <- scale(kronecker(rep(1,K) - ((1:K - K/2)^2)/K^2 + rnorm(K, 0, .1),
#'           1:K - K/2),
#'           center = FALSE, scale = TRUE)
#' tempb1 <- bs(x[, 1], K, intercept = TRUE)
#' tempb2 <- bs(x[, 2], K, intercept = TRUE)
#' tempW <- matrix(0, n, K^2)
#' for(i in 1:n) tempW[i,] <- kronecker(tempb1[i,], tempb2[i,])
#' def.x <- cbind(tempW%*%theta1, tempW%*%theta2)
#' fakeModel <- list(basis = list(B1 = bs(range(x[, 1]), K, intercept = TRUE),
#'                                B2 = bs(range(x[, 2]), K, intercept = TRUE)),
#'                   window = list(x = 0:1, y = 0:1),
#'                   x = x, def.x = def.x,
#'                   theta1 = theta1,
#'                   theta2 = theta2,
#'                   df1 = K,
#'                   df2 = K)
#'  plotGrid(fakeModel)
#'  plotGrid(fakeModel, margins = TRUE)
#'  plotGrid(fakeModel, persp = TRUE)
#'  # pdf("sim1.pdf"); plotGrid(fakeModel); dev.off()
#'  # pdf("sim2.pdf"); plotGrid(fakeModel, margins = TRUE); dev.off()
#'  # pdf("sim3.pdf"); plotGrid(fakeModel, persp = TRUE); dev.off()
#'
#' @author Guilherme Ludwig and Ronaldo Dias
#'
#' @references
#'
#'   To add.
#'
#' @seealso \code{\link{geoR::likfit}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
plotGrid <- function(model, nx = 20, ny = 20,
                     colorO = "Black", colorD = "Red",
                     plot = TRUE, margins = FALSE,
                     persp = FALSE, F1 = NULL, F2 = NULL) {
  theta1 <- model$theta1
  theta2 <- model$theta2
  def.x <- model$def.x
  rx <- model$window$x
  ry <- model$window$y
  if(is.null(model$basis)) {
    B <- list(B1 = bs(rx), B2 = bs(ry))
  } else {
    B <- model$basis
  }
  xyg <- model$x
  xgrid <- seq(rx[1], rx[2], length = nx)
  ygrid <- seq(ry[1], ry[2], length = ny)
  xygrid <- as.matrix(expand.grid(xgrid, ygrid))
  B1 <- predict(B$B1, xygrid[,1])
  B2 <- predict(B$B2, xygrid[,2])
  W <- matrix(0, nrow(xygrid), prod(sapply(B, dim)[2,]))
  for(i in 1:nrow(xygrid)) W[i,] <- kronecker(B1[i,], B2[i,])
  if(is.null(F1)) {
    f1 <- as.numeric(W%*%theta1)
  } else {
    f1 <- F1(xygrid[,1], xygrid[,2])
  }
  if(is.null(F2)) {
    f2 <- as.numeric(W%*%theta2)
  } else {
    f2 <- F2(xygrid[,1], xygrid[,2])
  }
  fgrid <- cbind(f1, f2)
  xl <- range(c(rx, f1))
  yl <- range(c(ry, f2))
  if(margins | persp) {
    layout(matrix(c(2,2,0,1,1,3,1,1,3), ncol = 3))
  }
  if(plot){
    plot(x, xlim = xl, ylim = yl, col = colorO,
         ylab = expression(y[2]), xlab = expression(y[1]))
    points(def.x, xlim = xl, ylim = yl, col = colorD)
    for (i in 1:ny) {
      lines(xygrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),], col = colorO)
      lines(fgrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),], col = colorD)
    }
    for (j in 1:nx) {
      lines(xygrid[seq((1+(j-1)),nx*ny,nx),], col = colorO)
      lines(fgrid[seq((1+(j-1)),nx*ny,nx),], col = colorD)
    }
    if(margins){
      image(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[nx:1,]),
            col = tim.colors(64),
            ylab = expression(x[2]), xlab = expression(x[1]))
      contour(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[nx:1,]),
              add = TRUE)
      image(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[ny:1,]),
            col = tim.colors(64),
            ylab = expression(x[2]), xlab = expression(x[1]))
      contour(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[ny:1,]),
              add = TRUE)
    }
    if(persp){
      persp(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[nx:1,]))
      persp(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[ny:1,]))
    }
  }
  invisible(fgrid)
}
