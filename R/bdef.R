#' Spatial Deformation via tensor product of B-splines
#'
#' This function finds a spatial deformation for Gaussian processes, based
#' on the penalized log-likelihood function. The deformation is obtained
#' via the tensor product of B-spline basis functions.
#'
#' @param cov.model A model for the spatial covariance function. See
#'                     \code{\link{RMmodel}} from the RandomFields package
#'                     for details. Defaults to Exponential model with
#'                     unknown variance and scale, plus a nugget effect. Right
#'                     now the model doesn't allow space-time covariance models.
#' @param x A n times 2 matrix of coordinates for the sampled data points. The
#'                     rows must be unique, corresponding to distinct sample
#'                     locations.
#' @param y A vector of length n*m with the values taken by the response
#'                     at the corresponding spatial locations in x, for
#'                     each time point t_1, ..., t_m.
#' @param tim A vector of equally spaced time positions, or NULL, in which
#'                     the temporal component of the model is ignored.
#'                     Defaults to NULL. Each time point must have the
#'                     same corresponding locations as rows in x.
#' @param df1 Degrees of freedom for the tensor product of splines along
#'                     the x[,1] coordinate. Defaults to 6.
#' @param df2 Degrees of freedom for the tensor product of splines along
#'                     the x[,2] coordinate. Defaults to 6.
#' @param window A list with two named entries, "x" and "y", each a length 2
#'                     numeric vector with the minimum and maximum values
#'                     taken for the respective coordinate in the sampled
#'                     region. Defaults to the range of the columns in
#'                     argument x.
#' @param maxit The maximum number of iterations between
#'                     optimization of the covariance parameters
#'                     and estimation of the spline coefficients.
#'                     Defaults to 2.
#' @param traceback Whether the function will return each step
#'                     of the optimization iterations, as a list.
#'                     For debugging purposes. Defaults to FALSE.
#' @param fullDes Whether the returned object will include a processed
#'                     copy of the dataset. Needed for the
#'                     \code{\link{plotGDdist}} function, but can be set
#'                     to FALSE if the user has no interest in GDdist.
#'                     Defaults to TRUE.
#' @param ... Additional arguments for RFfit.
#'
#' @export
#' @return \code{bdef} returns an object of \code{\link{class}} "bdef"
#' containing at least the following components
#'   \item{window}{The rectangular spatial domain in which the
#'                data was sampled}
#'   \item{basis}{list of B-spline basis functions, inherited from
#'                \code{\link{bs}}, with one set for each coordinate}
#'   \item{theta1}{B-spline coefficients for the coordinate x1}
#'   \item{theta2}{B-spline coefficients for the coordinate x2}
#'
#' @examples
#' # Example using artificially generated data
#' set.seed(1)
#' n <- 40
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x <- cbind(x1,x2)
#' F1 <- function(x1,x2) {
#'   x1
#' }
#' F2 <- function(x1,x2) {
#'   sqrt(x1*x2)
#' }
#' plotGrid(list(window = list(x = range(x1), y = range(x2)),
#'               x = x, def.x = cbind(F1(x1,x2), F2(x1,x2))),
#'          F1 = F1, F2 = F2, margins = TRUE)
#' covModel <- RMexp(var = 1, scale = 2) + RMnugget(var = 1)
#' data <- RFsimulate(covModel, x = F1(x1,x2), y = F2(x1,x2))
#' y <- as.numeric(unlist(data@data))
#' test1 <- bdef(x, y)
#' plotGrid(test1, margins = TRUE)
#' # \dontrun{
#' # test1 <- bdef(x, y, cov.model = RMmatern(nu = 2.5, var = NA, scale = NA))
#' # plotGrid(test2, margins = TRUE)
#' # }
#' # Time version
#' TIME <- 20
#' covModel2 <- RMexp(proj = "space", var = 1, scale = .25) *
#'   RMnugget(proj = "time", var = 1)
#' data2 <- RFsimulate(covModel2, x = F1(x1,x2), y = F2(x1,x2), T = 1:TIME)
#' y <- as.numeric(unlist(data2@data))
#' # \dontrun{
#' # matplot(1:TIME, t(matrix(y, ncol = TIME)), type = "b", lty = 1, pch = 1)
#' # }
#' covModel2m <- RMexp(var = NA, scale = NA) # Only space?
#' system.time({testT <- bdef(x, y, tim = 1:TIME, cov.model = covModel2m,
#'               maxit = 10, traceback = TRUE)})
#' plotGrid(testT, margins = TRUE)
#' plotGrid(list(basis = testT$basis,
#'               window = testT$window,
#'               x = testT$x,
#'               def.x = testT$x, # Fix here
#'               theta1 = testT$trace[[2]]$theta1,
#'               theta2 = testT$trace[[2]]$theta2,
#'               df1 = 6,
#'               df2 = 6), margins = TRUE)
#'               plotGrid(list(basis = testT$basis,
#'               window = testT$window,
#'               x = testT$x,
#'               def.x = testT$x, # Fix here
#'               theta1 = testT$trace[[4]]$theta1,
#'               theta2 = testT$trace[[4]]$theta2,
#'               df1 = 6,
#'               df2 = 6), margins = TRUE)
#'               plotGrid(list(basis = testT$basis,
#'               window = testT$window,
#'               x = testT$x,
#'               def.x = testT$x, # Fix here
#'               theta1 = testT$trace[[6]]$theta1,
#'               theta2 = testT$trace[[6]]$theta2,
#'               df1 = 6,
#'               df2 = 6), margins = TRUE)
#'
#' @author Guilherme Ludwig and Ronaldo Dias
#'
#' @references
#'
#'   To add.
#'
#' @seealso \code{\link{RandomFields}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
bdef <- function(x, y, tim = NULL,
                 cov.model = RMexp(var = NA, scale = NA) + RMnugget(var = NA),
                 df1 = 6, df2 = 6,
                 window = list(x = range(x[,1]), y = range(x[,2])),
                 maxit = 2, traceback = FALSE, fullDes = TRUE, ...) {

  RFoptions(printlevel=0, warn_normal_mode=FALSE)

  if(maxit < 0){
    maxit <- 1
  }
  # if(!is.null(tim) && length(tim) <= 3){
  #   # CHECK AMBIGUITY
  # }

  # Assuming total separability for now!
  #!# NEED TO SORT X!!
  n <- nrow(x)
  m <- length(tim)
  if(ncol(x)>2){
    tim <- unique(x[,3])
    x <- x[,-1*(3:ncol(x))]
    x <- unique(x)
  } else {
    x <- unique(x)
  }

  B1 <- bs(range(x[, 1]), df = df1, intercept = TRUE)
  B2 <- bs(range(x[, 2]), df = df2, intercept = TRUE)
  basis <- list(B1 = B1, B2 = B2)
  W <- matrix(0, n, df1*df2) # TODO debug unequal df
  for(i in 1:n) W[i,] <- kronecker(predict(basis$B1, x[, 1])[i, ],
                                   predict(basis$B2, x[, 2])[i, ])

  theta00 <- c(as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 1]))),
               as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 2]))))

  if(m == 0){
    model0 <- try(RFfit(model = cov.model, x = x[,1], y = x[,2], data = y, ...))
  } else {
    #!# model0 <- try(RFfit(model = cov.model, x = x[,1], y = x[,2], T = tim, data = y, ...))
    model0 <- try(RFfit(model = cov.model, x = x[,1], y = x[,2], data = matrix(y, nrow = n), ...))
  }

  if(m==0){
    theta0 <- auglag(par = theta00,
                     fn = likelihoodTarget, # gr = NULL
                     hin = jacobianConstraint, # hin.jac = NULL
                     control.outer = list(trace = FALSE,
                                          kkt2.check = FALSE),
                     DF1 = df1, DF2 = df2, B = basis,
                     M = model0, X = x, Y = y)$par
  } else {
    theta0 <- auglag(theta00,
                     fn = likelihoodTarget, # gr = NULL
                     hin = jacobianConstraint, # hin.jac = NULL
                     control.outer = list(trace = FALSE,
                                          kkt2.check = FALSE),
                     DF1 = df1, DF2 = df2, B = basis,
                     M = model0, X = x, Y = matrix(y, ncol = m))$par
  }

  # Traces the deformation map estimation
  if(traceback) {
    trace <- list()
    trace[["model_0"]] <- model0
    trace[["theta_0"]] <- list(theta1 = theta0[1:(df1*df2)],
                               theta2 = theta0[1:(df1*df2) + (df1*df2)])
  }

  f1 <- as.numeric(W%*%theta0[1:(df1*df2)])
  f2 <- as.numeric(W%*%theta0[1:(df1*df2) + (df1*df2)])

  if(m == 0){
    model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = y, ...))
    theta.new <- auglag(theta0,
                        fn = likelihoodTarget, # gr = NULL
                        hin = jacobianConstraint, # hin.jac = NULL
                        control.outer = list(trace = FALSE,
                                             kkt2.check = FALSE),
                        DF1 = df1, DF2 = df2, B = basis,
                        M = model0, X = x, Y = y)$par
  } else {
    # model1 <- try(RFfit(model = cov.model, x = f1, y = f2, T = tim, data = y, ...))
    model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = matrix(y, nrow = n), ...))
    theta.new <- auglag(theta0,
                        fn = likelihoodTarget, # gr = NULL
                        hin = jacobianConstraint, # hin.jac = NULL
                        control.outer = list(trace = FALSE,
                                             kkt2.check = FALSE),
                        DF1 = df1, DF2 = df2, B = basis,
                        M = model0, X = x, Y = matrix(y, ncol = m))$par
  }

  # Traces the deformation map estimation
  if(traceback) {
    trace[["model_1"]] <- model1
    trace[["theta_1"]] <- list(theta1 = theta.new[1:(df1*df2)],
                               theta2 = theta.new[1:(df1*df2) + (df1*df2)])
  }

  it <- 1
  condition <- TRUE
  while(it < maxit & condition){
    model0 <- model1
    theta0 <- theta.new
    f1 <- as.numeric(W%*%theta0[1:(df1*df2)])
    f2 <- as.numeric(W%*%theta0[1:(df1*df2) + (df1*df2)])
    if(m == 0){
      model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = y, ...))
      theta.new <- auglag(theta0,
                          fn = likelihoodTarget, # gr = NULL
                          hin = jacobianConstraint, # hin.jac = NULL
                          control.outer = list(trace = FALSE,
                                               kkt2.check = FALSE),
                          DF1 = df1, DF2 = df2, B = basis,
                          M = model0, X = x, Y = y)$par
    } else {
      # BUG HERE
      # model1 <- try(RFfit(model = cov.model, x = f1, y = f2, T = tim, data = y, ...))
      model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = matrix(y, nrow = n), ...))
      theta.new <- auglag(theta0,
                          fn = likelihoodTarget, # gr = NULL
                          hin = jacobianConstraint, # hin.jac = NULL
                          control.outer = list(trace = FALSE,
                                               kkt2.check = FALSE),
                          DF1 = df1, DF2 = df2, B = basis,
                          M = model0, X = x, Y = matrix(y, ncol = m))$par
    }
    condition <- all(max(abs(theta0 - theta.new)/abs(theta0), 0.0009) < 0.001)
    it <- it + 1
    # Traces the deformation map estimation
    if(traceback) {
      trace[[paste0("model_", it)]] <- model1
      trace[[paste0("theta_", it)]] <- list(theta1 = theta.new[1:(df1*df2)],
                                            theta2 = theta.new[1:(df1*df2) + (df1*df2)])
    }
  }

  ret <- list(window = window,
              basis = basis,
              x = x, def.x = cbind(f1, f2),
              theta1 = theta.new[1:(df1*df2)],
              theta2 = theta.new[1:(df1*df2) + (df1*df2)],
              df1 = df1,
              df2 = df2,
              model = model1)
  # Traces the deformation map estimation
  if(traceback){
    ret$trace <- trace
  }
  if(fullDes){
    if((!is.null(m)) && m > 1){
      ret$fullDes <- list(x = cbind(rep(x[,1], m), rep(x[,2], m)),
                          def.x = cbind(rep(f1, m), rep(f2, m)),
                          tim = tim,
                          y = y)
    } else {
      ret$fullDes <- list(x = x,
                          def.x = cbind(f1,f2),
                          tim = tim,
                          y = y)
    }
  }

  class(ret) <- "bdef"

  RFoptions(printlevel=1, warn_normal_mode=TRUE)

  return(ret)

}
