#' Plot deformation traceback as a function of iterations
#'
#' This function is a debugging tool.
#'
#' @export
#' @examples
#' # Example using artificially generated data
#' set.seed(1)
#' n <- 50
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x <- cbind(x1,x2)
#' F1 <- function(x1,x2) {
#'   x1 - .2*sqrt(x1*x2)
#' }
#' F2 <- function(x1,x2) {
#'   .8*x2 + .4*sqrt(x1*x2) - .6*sqrt((1-x2)*x1)
#' }
#' TIME <- 20
#' covModel2 <- RMexp(proj = "space", var = .8, scale = 1) + RMnugget(proj = "time", var = .2)
#' data2 <- RFsimulate(covModel2, x = F1(x1,x2), y = F2(x1,x2), T = 1:TIME)
#' y <- as.numeric(unlist(data2@data))
#' covModel2m <- RMexp(var = NA, scale = NA) + RMnugget(var=NA) # Only space, independent in time
#' system.time({testT <- bdef(x, y, tim = 1:TIME, cov.model = covModel2m,
#'               maxit = 10, traceback = TRUE)})
#' plotGrid(testT, margins = TRUE)
#' plotTraceback(testT, margins = TRUE)
plotTraceback <- function(){

  plotGrid(list(basis = testT$basis,
                window = testT$window,
                x = testT$x,
                def.x = testT$x, # Fix here
                theta1 = testT$trace[[2]]$theta1,
                theta2 = testT$trace[[2]]$theta2,
                df1 = 6,
                df2 = 6), margins = TRUE)
  plotGrid(list(basis = testT$basis,
                window = testT$window,
                x = testT$x,
                def.x = testT$x, # Fix here
                theta1 = testT$trace[[4]]$theta1,
                theta2 = testT$trace[[4]]$theta2,
                df1 = 6,
                df2 = 6), margins = TRUE)
  plotGrid(list(basis = testT$basis,
                window = testT$window,
                x = testT$x,
                def.x = testT$x, # Fix here
                theta1 = testT$trace[[6]]$theta1,
                theta2 = testT$trace[[6]]$theta2,
                df1 = 6,
                df2 = 6), margins = TRUE)

}
