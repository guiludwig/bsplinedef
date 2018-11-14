#' Find Tuning parameters for B-spline spatial deformation
#' 
#' This auxiliary function helps finding tuning parameters (lambda, zeta) for
#' the B-spline tensor product spatial deformation. Execution might be slow.
#'
#' @param model an object from the bdef function.
#' @param grid.lambda Grid of values for the penalty parameter lambda.
#' @param grid.zeta Grid of values for the penalty parameter zeta.
#' @param verbose Verbose execution informs user of the LOOCV progress. Defaults to FALSE.
#'
#' @export
#' 
#' @return a data frame with the selected grid, and the leave-one-out cross validation error
#' 
#' @examples
#' # Example using artificially generated data
#' \dontrun{
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
#' covModel <- RMexp(var = 1, scale = .25, proj = "space") + RMnugget(var = 1) # Independent in time
#' data <- RFsimulate(covModel, x = F1(x[,1],x[,2]), y = F2(x[,1],x[,2]), 
#'                    T = seq(from = 1, by = 1, len = TIME)) # order ~ expand.grid(x, y, T)
#' y <- as.numeric(unlist(data@data))
#' # Model for spatial dependence, time is assumed independent
#' covModelM <- RMexp(var = NA, scale = NA) + RMnugget(var = NA)
#' # Calculates deformation, profle likelihood up to maxit times
#' test.def <- bdef(x, y, tim = 1:TIME, cov.model = covModelM, maxit = 10)
#' tuningChoice <- findTuning(test.def, verbose = TRUE)
#' require(ggplot2)
#' require(dplyr)
#' tuningChoice %>%
#' mutate(zeta = factor(zeta)) %>%
#' mutate(lambda = log(lambda)) %>%
#' ggplot(aes(x = lambda, y = LOOCV, color = zeta)) + 
#' geom_point() +
#' geom_line(aes(group = zeta))
#' }
findTuning <- function(model, 
                       grid.lambda = 10^(-5:5),
                       grid.zeta = min(apply(model$x, 2, function(x) diff(range(x))))*(1:10)/20,
                       verbose = FALSE){
  stopifnot(class(model) == "bdef")
  results <- expand.grid(lambda = grid.lambda, zeta = grid.zeta)
  results$LOOCV <- numeric(nrow(results))
  for(param in 1:nrow(results)){
    for(i in 1:nrow(model$x)){
      rmData <- seq(i, length(model$y), nrow(model$x))
      model.new <- bdef(x = model$x[-i,], y = model$y[-rmData], tim = model$tim,
                        cov.model = model$cov.model,
                        df1 = model$df1, df2 = model$df2,
                        lambda = results$lambda[param], 
                        zeta1 = results$zeta[param], zeta2 = results$zeta[param],
                        window = model$window)
      results$LOOCV[param] <- results$LOOCV[param] + sum((as.numeric(predict(model.new, model$x[i, , drop = FALSE])$krige) - 
                                                            model$y[rmData])^2)/length(rmData)
      if(verbose) cat("Iteration ", param, " of ", nrow(results), ", step ", i, " of ", nrow(model$x),".\n")
    }
  }
  return(results)
}