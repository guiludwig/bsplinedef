#' Plot deformation traceback as a function of iterations
#'
#' This function is a debugging tool.
#'
#' @export
#' @examples
#' # Example using artificially generated data
#' # Example using artificially generated data
#' set.seed(1)
#' m <- 10
#' x1 <- (0:m)/m
#' x2 <- (0:m)/m
#' x <- as.matrix(expand.grid(x1,x2))
#' n <- nrow(x)
#' F1 <- function(x1,x2, a = 2.5, b = 1.0) {
#' x <- x1 - 0.5; y <- x2 - 0.5
#' angle <- a*exp(-(x*x+y*y)/(b*b)) + 3*pi/2
#' return(cos(angle)*x + sin(angle)*y + 0.5)
#' }
#' F2 <- function(x1,x2, a = 2.5, b = 1.0) {
#' x <- x1 - 0.5; y <- x2 - 0.5
#' angle <- a*exp(-(x*x+y*y)/(b*b)) + 3*pi/2
#' return(-sin(angle)*x + cos(angle)*y + 0.5)
#' }
#' TIME <- 20
#' covModel <- RMexp(var = 1, scale = .25) + RMnugget(var = 1) # Independent in time
#' data <- RFsimulate(covModel, x = F1(x[,1],x[,2]), y = F2(x[,1],x[,2]), T = 1:TIME)
#' y <- as.numeric(unlist(data@data))
#' covModelM <- RMexp(var = NA, scale = NA) + RMnugget(var = NA)
#' test.def <- bdef(x, y, tim = 1:TIME, cov.model = covModelM, maxit = 10)
#' # Estimated deformation
#' plotGrid(test.def)
#' plotTraceback(test.def, margins = TRUE)
plotTraceback <- function(model, sleep = 2, ...){

  if(is.null(model$trace)){
    stop("Enable traceback when fitting the model, see function bdef().")
  }
  dev.hold()
  plotGrid(list(basis = model$basis,
                window = model$window,
                x = model$x,
                def.x = model$x, # Fix here
                theta1 = model$trace[["theta_1"]]$theta1,
                theta2 = model$trace[["theta_1"]]$theta2,
                df1 = model$df1,
                df2 = model$df1), ...)
  dev.flush()
  p <- 2
  while(!is.null(model$trace[[paste0("theta_",p)]])){
    dev.hold()
    plotGrid(list(basis = model$basis,
                  window = model$window,
                  x = model$x,
                  def.x = model$x, # Fix here
                  theta1 = model$trace[["theta_1"]]$theta1,
                  theta2 = model$trace[["theta_1"]]$theta2,
                  df1 = model$df1,
                  df2 = model$df1), ...)
    dev.flush()
  }

}