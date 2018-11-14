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
#'                     shown along the corresponding axes. Defaults to TRUE.
#' @param F1 Allows the user to supply a function that generates the deformed
#'                     coordinate y1. Defaults to NULL, in which case the
#'                     tensor product of B-splines is used instead.
#' @param F2 Allows the user to supply a function that generates the deformed
#'                     coordinate y2. Defaults to NULL, in which case the
#'                     tensor product of B-splines is used instead.
#' @param chull Whether the deformation is displayed inside of the convex hull
#'                     of the sampled sites, or a square window containing the sampled
#'                     sites. Defaults to FALSE, which corresponds to the square
#'                     window.
#' @param ... Extra parameters to be passed to perspective plot, in particular the
#'                     azymutal direction and colatitude. See \code{\link{persp}} for
#'                     a description.
#'
#' @export
#' @return \code{plotGrid} invisibly returns the values of the deformation on
#'                     the grid.
#'
#' @examples
#' # Example using artificially generated data
#' set.seed(1)
#' m <- 10
#' x1 <- (0:m)/m
#' x2 <- (0:m)/m
#' x <- as.matrix(expand.grid(x1,x2))
#' n <- nrow(x)
#' K <- 4
#' tempb1 <- bs(x[, 1], K, intercept = TRUE)
#' tempb2 <- bs(x[, 2], K, intercept = TRUE)
#' tempW <- matrix(0, n, K^2)
#' for(i in 1:n) tempW[i,] <- kronecker(tempb2[i,], tempb1[i,])
#' theta1 <- solve(crossprod(tempW), crossprod(tempW, x[,1] + .2*x[,2]))
#' theta2 <- solve(crossprod(tempW), crossprod(tempW, x[,2] - .2*sqrt(x[,1]*x[,2])))
#' def.x <- cbind(tempW%*%theta1, tempW%*%theta2)
#' fakeModel <- list(basis = list(B1 = bs(range(x[, 1]), K, intercept = TRUE),
#'                                B2 = bs(range(x[, 2]), K, intercept = TRUE)),
#'                   window = list(x = 0:1, y = 0:1),
#'                   x = x, def.x = def.x,
#'                   theta1 = theta1,
#'                   theta2 = theta2,
#'                   df1 = K,
#'                   df2 = K)
#'  plotGrid(fakeModel, margins = FALSE)
#'  plotGrid(fakeModel)
#'  plotGrid(fakeModel, persp = TRUE, theta = 15)
#' # Example with truth function
#' set.seed(1)
#' m <- 10
#' x1 <- (0:m)/m
#' x2 <- (0:m)/m
#' x <- as.matrix(expand.grid(x1,x2))
#' n <- nrow(x)
#' F1 <- function(x1,x2, a = 2.5, b = 1.0) {
#'   x <- x1 - 0.5; y <- x2 - 0.5
#' angle <- a*exp(-(x*x+y*y)/(b*b)) + 3*pi/2
#' return(cos(angle)*x + sin(angle)*y + 0.5)
#' }
#' F2 <- function(x1,x2, a = 2.5, b = 1.0) {
#'   x <- x1 - 0.5; y <- x2 - 0.5
#' angle <- a*exp(-(x*x+y*y)/(b*b)) + 3*pi/2
#' return(y + x^2)
#' }
#' plotGrid(list(window = list(x = range(x1), y = range(x2)),
#' x = x, def.x = cbind(F1(x[,1],x[,2]), 
#'                                   F2(x[,1],x[,2]))),
#'         F1 = F1, F2 = F2)
#'
#' @author Guilherme Ludwig and Ronaldo Dias
#'
#' @references
#'
#'   To add.
#'
#' @seealso \code{\link{bdef}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
plotGrid <- function(model, nx = 20, ny = 20,
                     colorO = "Black", colorD = "Red",
                     plot = TRUE, margins = TRUE,
                     persp = FALSE, F1 = NULL, F2 = NULL,
                     chull = FALSE, offset = c(0,0), ...) {
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
  chx <- chull(xyg) # Convex Hull, clockwise locations of points
  centerHull <- mean(xyg[chx,]) # Average = baricenter of hull?
  chx <- c(chx, chx[1]) # Loops around the hull, if needed
  xgrid <- seq(rx[1], rx[2], length = nx)
  ygrid <- seq(ry[1], ry[2], length = ny)
  xygrid <- as.matrix(expand.grid(xgrid, ygrid))
  B1 <- predict(B$B1, xygrid[,1])
  B2 <- predict(B$B2, xygrid[,2])
  W <- matrix(0, nrow(xygrid), prod(sapply(B, dim)[2,]))
  for(i in 1:nrow(xygrid)) W[i,] <- kronecker(B2[i,], B1[i,])
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
    points(sweep(def.x, 2, offset[1:2]), xlim = xl, ylim = yl, col = colorD)
    if(chull){
      for (i in 1:ny) {
        xygrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),]
        tempN <- nrow(temp)
        condition <- rep(TRUE, tempN)
        for(j in 1:tempN){
          if(condition[j]){
            lines(xygrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),], col = colorO)
            lines(fgrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),], col = colorD)
          }
        }
      }
      for (j in 1:nx) {
        temp <- xygrid[seq((1+(j-1)),nx*ny,nx),]
        tempN <- nrow(temp)
        condition <- rep(TRUE, tempN) # is
        for(i in 1:tempN){
          if(condition[i]){
            lines(xygrid[seq((1+(j-1)),nx*ny,nx),], col = colorO)
            lines(fgrid[seq((1+(j-1)),nx*ny,nx),], col = colorD)
          }
        }
      }
    } else {
      for (i in 1:ny) {
        lines(xygrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),], col = colorO)
        lines(fgrid[(1+(i-1)*nx):((1+(i-1)*nx)+nx-1),]-offset[1], col = colorD)
      }
      for (j in 1:nx) {
        lines(xygrid[seq((1+(j-1)),nx*ny,nx),], col = colorO)
        lines(fgrid[seq((1+(j-1)),nx*ny,nx),]-offset[2], col = colorD)
      }
    }
    if(margins){
      image(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[1:nx,]),
            col = tim.colors(64),
            ylab = expression(x[2]), xlab = expression(x[1]))
      contour(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[1:nx,]),
              add = TRUE)
      image(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[1:ny,]),
            col = tim.colors(64),
            ylab = expression(x[2]), xlab = expression(x[1]))
      contour(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[1:ny,]),
              add = TRUE)
    }
    if(persp){
      # Is this dimension flip correct? Also persp is pretty bad to visualize
      persp(list(x = xgrid, y = ygrid, z = matrix(f2, nx, ny)[1:nx,]),
            xlab = expression(x[1]), ylab = expression(x[2]), zlab = expression(y[1]),
            ...)
      persp(list(x = ygrid, y = xgrid, z = matrix(f1, nx, ny)[1:ny,]),
            xlab = expression(x[1]), ylab = expression(x[2]), zlab = expression(y[1]),
            ...)
    }
  }
  invisible(fgrid)
}
