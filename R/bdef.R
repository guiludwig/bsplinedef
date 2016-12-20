#' Spatial Deformation via B-spline
#'
#' This function finds a spatial deformation for Gaussian processes, based
#' on the penalized log-likelihood function. The deformation is obtained
#' via the tensor product of B-spline basis functions.
#'
#' @param cov.model A model for the covariance function. See \code{\link{RMmodel}}
#'                     on the RandomFields package for details. Defaults to
#'                     Mat\'{e}rn model with parameters nu = 2.5 and
#'                     variance and scale to be estimated.
#' @param x A n times 2 matrix of coordinates for the sampled data points.
#' @param y A vector of length n with the values taken by the response
#'                     at the corresponding spatial locations in x.
#' @param def.x A n times 2 matrix of deformed coordinates,
#' @param df1 Square root of the total degrees of freedom for the tensor
#'                     product of B-spline basis for the first coordinate.
#'                     Defaults to 8.
#' @param df2 Square root of the total degrees of freedom for the tensor
#'                     product of B-spline basis for the second coordinate.
#'                     Defaults to 8.
#' @param window A list with two named entries, "x" and "y", each a length 2
#'                     numeric vector with the minimum and maximum values
#'                     taken for the respective coordinate in the sampled
#'                     region. Defaults to the range of the columns in
#'                     argument x.
#'  @param lambda Penalization parameter. If set to NA (default), this
#'                     parameter will be estimated. Otherwise, it will use
#'                     the user supplied value. The parameter is in a log
#'                     scale, thus ranging from -Inf to +Inf.
#'  @param maxit The maximum number of iterations between
#'                     optimization of the covariance parameters
#'                     and estimation of the spline coefficients.
#'                     Defaults to 2.
#'  @param traceback Whether the function will return each step
#'                     of the optimization iterations, as a list.
#'                     For debugging purposes. Defaults to FALSE.
#'  @param ... Additional arguments for RFfit.
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
#'   x1 # (x1 - 2*x2)^2/sqrt(x1^2 + 4*x2^2)
#' }
#' F2 <- function(x1,x2) {
#'   sqrt(x1*x2) # (3*x1 + x2)^2/sqrt(9*x1^2 + x2^2)
#' }
#' plotGrid(list(window = list(x = range(x1), y = range(x2)),
#'               x = x, def.x = cbind(F1(x1,x2), F2(x1,x2))),
#'          F1 = F1, F2 = F2, margins = TRUE)
#' covModel <- RMmatern(nu = 2.5, var = 1, scale = 2) # No nugget
#' data <- RFsimulate(covModel, x = F1(x1,x2), y = F2(x1,x2))
#' y <- as.numeric(unlist(data@data))
#' test1 <- bdef(x, y) # lambda estimated
#' test2 <- bdef(x, y, lambda = -100) # lambda provided by user
#' test3 <- bdef(x, y, lambda = 100) # ditto
#' test4 <- bdef(x, y, cov.model = RMmatern(nu = 2.5, var = NA, scale = NA)) # No nugget effect
#' plotGrid(test1, margins = TRUE)
#' plotGrid(test2, margins = TRUE)
#' plotGrid(test3, margins = TRUE)
#' plotGrid(test4, margins = TRUE)
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
bdef <- function(x, y,
                 cov.model = RMmatern(nu = 2.5, var = NA, scale = NA) + RMnugget(var = NA),
                 df1 = 6, df2 = 6,
                 window = list(x = range(x[,1]), y = range(x[,2])),
                 lambda = NA, maxit = 2,
                 traceback = FALSE, ...) {

  RFoptions(printlevel=0, warn_normal_mode=FALSE)

  if(maxit < 0){
    maxit <- 1
  }

  B1 <- bs(range(x[, 1]), df = df1, intercept = TRUE)
  B2 <- bs(range(x[, 2]), df = df2, intercept = TRUE)
  basis <- list(B1 = B1, B2 = B2)
  W <- matrix(0, n, df1*df2) # TODO debug unequal df
  for(i in 1:n) W[i,] <- kronecker(predict(basis$B1, x[, 1])[i, ],
                                   predict(basis$B2, x[, 2])[i, ])
  theta00 <- c(as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 1]))),
               as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 2]))))

  model0 <- try(RFfit(model = cov.model, x = x[,1], y = x[,2], data = y))

  # bsplinedef:::likelihoodTarget(c(seq(1, 2*df1*df2), 1))
  if(is.na(lambda)){
    # theta0 <- optim(c(seq(1, 2*df1*df2), 1), likelihoodTarget)$par
    theta0 <- optim(c(theta00, 0), likelihoodTarget,
                    DF1 = df1, DF2 = df2, B = basis,
                    M = model0, X = x, Y = y, L = lambda)$par
  } else {
    theta0 <- optim(theta00, likelihoodTarget, # seq(1, 2*df1*df2)
                    DF1 = df1, DF2 = df2, B = basis,
                    M = model0, X = x, Y = y, L = lambda)$par
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
  if(is.na(lambda)){
    hat.lambda = theta0[2*df1*df2+1]
  } else {
    hat.lambda = lambda
  }

  model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = y))
  theta.new <- optim(theta0, likelihoodTarget,
                     DF1 = df1, DF2 = df2, B = basis,
                     M = model1, X = x, Y = y, L = lambda)$par

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
    model1 <- try(RFfit(model = cov.model, x = f1, y = f2, data = y))
    theta.new <- optim(theta0, likelihoodTarget,
                       DF1 = df1, DF2 = df2, B = basis,
                       M = model1, X = x, Y = y, L = lambda)$par
    condition <- all(max(abs(theta0 - theta.new)/theta0, 0.0009) < 0.001)
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
              lambda = hat.lambda,
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
  class(ret) <- "bdef"

  RFoptions(printlevel=1, warn_normal_mode=TRUE)

  return(ret)

}
