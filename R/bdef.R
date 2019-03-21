#' Spatial Deformation via tensor product of B-splines
#'
#' This function finds a spatial deformation for Gaussian processes, based
#' on the penalized log-likelihood function. The deformation is obtained
#' via the tensor product of B-spline basis functions.
#'
#' @param cov.model A model for the spatial covariance function. See
#'                     \code{\link{RMmodel}} from the RandomFields package
#'                     for details. Defaults to Exponential model with
#'                     unknown variance and scale, without a nugget effect. Right
#'                     now the model doesn't allow space-time covariance models.
#'                     TODO: default RMexp(var = NA, scale = NA) + RMnugget(var = NA)
#' @param x A n times 2 matrix of coordinates for the sampled data points. The
#'                     rows must be unique, corresponding to distinct sample
#'                     locations. Alternatively, a n times 3 matrix with the third
#'                     column corresponding to time. The times must be in a regular
#'                     grid and exist in all spatial locations (TODO generalize).
#' @param y A vector of length n*m with the values taken by the response
#'                     at the corresponding spatial locations in x, for
#'                     each time point t_1, ..., t_m.
#' @param tim A vector of equally spaced time positions, or NULL, in which
#'                     the temporal component of the model is ignored.
#'                     Defaults to NULL. Each time point must have the
#'                     same corresponding locations as rows in x.
#' @param df1 Degrees of freedom for the tensor product of splines along
#'                     the x[,1] coordinate. Defaults to 4.
#' @param df2 Degrees of freedom for the tensor product of splines along
#'                     the x[,2] coordinate. Defaults to 4.
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
#' @param type Method to fit the deformation model, a choice of "jacobian" (which
#'                     includes a regularization term for the positive jacobian, slow),
#'                     "penalized" (quadratic penalty, default) and "none" (no penalty).
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
#' # No deformation reference, entries are independent in time
#' test.nondef <- RFfit(covModelM, x = x[,1], y = x[,2], 
#'                      # T = seq(from = 1, by = 1, len = TIME), # slow?
#'                      # data = y)
#'                      data = matrix(y, ncol = TIME))
#' # Calculates deformation, profle likelihood up to maxit times
#' test.def <- bdef(x, y, tim = 1:TIME, cov.model = covModelM, maxit = 10)
#' # Estimated deformation
#' plotGrid(test.def)
#' # Comparison of Variograms
#' plot(test.nondef, ylim = c(0,2), xlim = c(0,0.7),
#'      model = list(`true model` = RMexp(var = 1, scale = .25) + RMnugget(var = 1)))
#' plot(test.def$model, ylim = c(0,2), xlim = c(0,0.7),
#'      model = list(`true model` = RMexp(var = 1, scale = .25) + RMnugget(var = 1)))
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
                 type = c("penalized", "jacobian", "none"),
                 target = c("likelihood", "variogram"),
                 df1 = 6, df2 = 6, lambda = .5, 
                 zeta1 = .5, zeta2 = .5,
                 window = list(x = range(x[,1]), y = range(x[,2])),
                 maxit = 2, traceback = TRUE, 
                 fullDes = TRUE, 
                 ...) {
  
  RFoptions(printlevel = 0, warn_normal_mode = FALSE)
  
  type <- match.arg(type)
  target <- match.arg(target)
  
  if(maxit < 0){
    maxit <- 1
  }
  if(!is.matrix(x)){
    message("x coerced to be a matrix, please check input.")
    x <- as.matrix(x)
  }
  
  # Assuming total separability for now!
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
  dB1 <- splineDesign(knots = c(rep(attr(B1, "Boundary.knots")[1],
                                    attr(B1, "degree") + 1),
                                attr(B1, "knots"),
                                rep(attr(B1, "Boundary.knots")[2],
                                    attr(B1, "degree") + 1)),
                      x = x[, 1, drop = FALSE], outer.ok = TRUE,
                      derivs = 1)
  dB2 <- splineDesign(knots = c(rep(attr(B2, "Boundary.knots")[1],
                                    attr(B2, "degree") + 1),
                                attr(B2, "knots"),
                                rep(attr(B2, "Boundary.knots")[2],
                                    attr(B2, "degree") + 1)),
                      x = x[, 2, drop = FALSE], outer.ok = TRUE,
                      derivs = 1)
  B1 <- predict(B1, x[ , 1, drop = FALSE])
  B2 <- predict(B2, x[ , 2, drop = FALSE])
  
  W <- matrix(0, n, df1*df2)
  for(i in 1:n) W[i,] <- kronecker(B2[i, ], B1[i, ])
  
  theta00 <- c(as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 1]))),
               as.numeric(solve(crossprod(W) + .001*diag(df1*df2)/max(W), crossprod(W, x[, 2]))))
  
  suppressMessages({
    model0 <- try(RFfit(model = cov.model, x = x[,1], y = x[,2], 
                        data = matrix(y, nrow = n), ...),
                  silent = TRUE)
  })
  
  if(target == "likelihood"){
    theta0 <- switch(type, 
                     penalized = optim(theta00,
                                       fn = likelihoodTargetPen, # gr = dLikelihoodTarget,
                                       DF1 = df1, DF2 = df2,
                                       M = model0, X = x, w = W,
                                       Y = matrix(y, nrow = n),
                                       LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par,
                     jacobian = auglag(theta00,
                                       fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                       hin = jacobianConstraint, # hin.jac = dJacobianConstraint,
                                       control.outer = list(trace = FALSE,
                                                            kkt2.check = FALSE),
                                       DF1 = df1, DF2 = df2,
                                       b1 = B1, b2 = B2,
                                       db1 = dB1, db2 = dB2,
                                       M = model0, X = x, w = W,
                                       Y = matrix(y, ncol = m))$par,
                     none = optim(theta00,
                                  fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                  DF1 = df1, DF2 = df2,
                                  M = model0, X = x, w = W,
                                  Y = matrix(y, nrow = n))$par)
  } else {
    #!#  nMDS: Sampson & Guttorp p. 110
    Sigma <- RFcovmatrix(cov.model, x = x[,1], y = x[,2])
    SV <- kronecker(matrix(1, nrow = nrow(Sigma)), 
                    matrix(diag(Sigma), ncol = ncol(Sigma))) + 
      kronecker(matrix(1, ncol = ncol(Sigma)), 
                matrix(diag(Sigma), nrow = nrow(Sigma))) +
      -2*Sigma # Sample Variog
    hatf0 <- cmdscale(SV) # TODO: Nonmetric, instead of metric
    theta0 <- optim(theta00,
                     fn = mdsTargetPen, 
                     DF1 = df1, DF2 = df2,
                     M = model0, X = x, w = W,
                     Y = hatf0,
                     LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par
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
  
  suppressMessages({
    model1 <- try(RFfit(model = cov.model, x = f1, y = f2, # T = tim, 
                        data = matrix(y, nrow = n), ...), silent = TRUE)
  })
  
  if(target == "likelihood"){
    theta.new <- switch(type, 
                        penalized = optim(theta0,
                                          fn = likelihoodTargetPen, # gr = dLikelihoodTarget,
                                          DF1 = df1, DF2 = df2,
                                          M = model0, X = x, w = W,
                                          Y = matrix(y, nrow = n),
                                          LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par,
                        jacobian = auglag(theta0,
                                          fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                          hin = jacobianConstraint, # hin.jac = dJacobianConstraint,
                                          control.outer = list(trace = FALSE,
                                                               kkt2.check = FALSE),
                                          DF1 = df1, DF2 = df2,
                                          b1 = B1, b2 = B2,
                                          db1 = dB1, db2 = dB2,
                                          M = model0, X = x, w = W,
                                          Y = matrix(y, ncol = m))$par,
                        none = optim(theta0,
                                     fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                     DF1 = df1, DF2 = df2,
                                     M = model0, X = x, w = W,
                                     Y = matrix(y, nrow = n))$par)
  } else {
    #!#  nMDS: Sampson & Guttorp p. 110
    Sigma <- RFcovmatrix(cov.model, x = f1, y = f2)
    SV <- kronecker(matrix(1, nrow = nrow(Sigma)), 
                    matrix(diag(Sigma), ncol = ncol(Sigma))) + 
      kronecker(matrix(1, ncol = ncol(Sigma)), 
                matrix(diag(Sigma), nrow = nrow(Sigma))) +
      -2*Sigma # Sample Variog
    hatf0 <- cmdscale(SV) # TODO: Nonmetric, instead of metric
    theta.new <- optim(theta0,
                    fn = mdsTargetPen, 
                    DF1 = df1, DF2 = df2,
                    M = model0, X = x, w = W,
                    Y = hatf0,
                    LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par
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
    
    suppressMessages({
      model1 <- try(RFfit(model = cov.model, x = f1, y = f2, # T = tim, 
                          data = matrix(y, nrow = n), ...), silent = TRUE)
    })
    if(target == "likelihood"){
      theta.new <- switch(type, 
                          penalized = optim(theta0,
                                            fn = likelihoodTargetPen, # gr = dLikelihoodTarget,
                                            DF1 = df1, DF2 = df2,
                                            M = model0, X = x, w = W,
                                            Y = matrix(y, nrow = n),
                                            LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par,
                          jacobian = auglag(theta0,
                                            fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                            hin = jacobianConstraint, # hin.jac = dJacobianConstraint,
                                            control.outer = list(trace = FALSE,
                                                                 kkt2.check = FALSE),
                                            DF1 = df1, DF2 = df2,
                                            b1 = B1, b2 = B2,
                                            db1 = dB1, db2 = dB2,
                                            M = model0, X = x, w = W,
                                            Y = matrix(y, ncol = m))$par,
                          none = optim(theta0,
                                       fn = likelihoodTarget, # gr = dLikelihoodTarget,
                                       DF1 = df1, DF2 = df2,
                                       M = model0, X = x, w = W,
                                       Y = matrix(y, nrow = n))$par)
    } else {
      #!#  nMDS: Sampson & Guttorp p. 110
      Sigma <- RFcovmatrix(cov.model, x = f1, y = f2)
      SV <- kronecker(matrix(1, nrow = nrow(Sigma)), 
                      matrix(diag(Sigma), ncol = ncol(Sigma))) + 
        kronecker(matrix(1, ncol = ncol(Sigma)), 
                  matrix(diag(Sigma), nrow = nrow(Sigma))) +
        -2*Sigma # Sample Variog
      hatf0 <- cmdscale(SV) # TODO: Nonmetric, instead of metric
      theta.new <- optim(theta0,
                         fn = mdsTargetPen, 
                         DF1 = df1, DF2 = df2,
                         M = model0, X = x, w = W,
                         Y = hatf0,
                         LAMBDA = lambda, z1 = zeta1, z2 = zeta2)$par
    }
    tempCon <- abs(theta0 - theta.new)/abs(theta0)
    tempCon <- ifelse(is.nan(tempCon), TRUE, tempCon < 1e-6) # ignores theta0 = 0 entries, keep opt
    condition <- all(tempCon)
    it <- it + 1
    # Traces the deformation map estimation
    if(traceback) {
      trace[[paste0("model_", it)]] <- model1
      trace[[paste0("theta_", it)]] <- list(theta1 = theta.new[1:(df1*df2)],
                                            theta2 = theta.new[1:(df1*df2) + (df1*df2)])
    }
  }
  
  f1 <- as.numeric(W%*%theta.new[1:(df1*df2)])
  f2 <- as.numeric(W%*%theta.new[1:(df1*df2) + (df1*df2)])
  
  ret <- list(window = window,
              basis = basis,
              x = x, def.x = cbind(f1, f2), y = y, tim = tim,
              theta1 = theta.new[1:(df1*df2)],
              theta2 = theta.new[1:(df1*df2) + (df1*df2)],
              df1 = df1,
              df2 = df2,
              model = model1, 
              cov.model = cov.model,
              type = type,
              target = target)
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
  
  RFoptions(printlevel = 1, warn_normal_mode = TRUE)
  
  return(ret)
  
}