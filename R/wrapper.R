#' MAnifold-constrained Gaussian process Inference (MAGI)
#'
#' @description Core function of the MAGI method for inferring the parameters and trajectories of dynamic systems governed by ordinary differential equations.
#' 
#' @param y data matrix of observations
#' @param odeModel list of ODE functions and inputs. See details.
#' @param tvec vector of discretization time points corresponding to rows of \code{y}.  If missing, \code{MagiSolver} will use the column named `time` in \code{y}.
#' @param control list of control variables, which may include `sigma`, `phi`, `theta`, `xInit`, `mu`, `dotmu`, `priorTemperature`, `niterHmc`,
#'  `nstepsHmc`, `burninRatio`, `stepSizeFactor`, `bandSize`, `useFixedSigma`, `kerneltype`, `skipMissingComponentOptimization`, `positiveSystem`, `verbose`.  See details.
#' 
#' @return 
#' \code{MagiSolver} returns an object of class \code{\link{magioutput}} which contains the following elements:
#' \describe{
#' \item{\code{theta}}{matrix of MCMC samples for the system parameters \eqn{\theta}, after burn-in.}
#' \item{\code{xsampled}}{array of MCMC samples for the system trajectories at each discretization time point, after burn-in.}
#' \item{\code{sigma}}{matrix of MCMC samples for the observation noise SDs \eqn{\sigma}, after burn-in.}
#' \item{\code{phi}}{matrix of estimated GP hyper-parameters, one column for each system component.}
#' \item{\code{lp}}{vector of log-posterior values at each MCMC iteration, after burn-in.}
#' \item{\code{y, tvec, odeModel}}{from the inputs to \code{MagiSolver}.}
#' }
#' 
#' @details 
#' The data matrix \code{y} has a column for each system component, and optionally a column `time` with the discretization time points. If the column `time` is not provided in \code{y}, a vector of time points must be provided via the \code{tvec} argument. The rows of \code{y} correspond to the discretization set \eqn{I} at which the GP is constrained to the derivatives of the ODE system. To set the desired discretization level for inference, use \code{\link{setDiscretization}} to prepare the data matrix for input into \code{MagiSolver}. Missing observations are indicated with \code{NA} or \code{NaN}.
#' 
#' The list \code{odeModel} is used for specification of the ODE system and its parameters. It must include five elements:
#' \describe{
#' \item{\code{fOde}}{function that computes the ODEs, specified with the form \code{f(theta, x, tvec)}. \code{fOde} should return a matrix where columns correspond to the system components of \code{x}, see examples.}
#' \item{\code{fOdeDx}}{function that computes the gradients of the ODEs with respect to the system components. \code{fOdeDx} should return a 3-D array, where the slice \code{[, i, j]} is the partial derivative of the ODE for the j-th system component with respect to the i-th system component, see examples.}
#' \item{\code{fOdeDtheta}}{function that computes the gradients of the ODEs with respect to the parameters \eqn{\theta}. \code{fOdeDtheta} should return a 3-D array, where the slice \code{[, i, j]} is the partial derivative of the ODE for the j-th system component with respect to the i-th parameter in \eqn{\theta}, see examples.}
#' \item{\code{thetaLowerBound}}{a vector indicating the lower bounds of each parameter in \eqn{\theta}.}
#' \item{\code{thetaUpperBound}}{a vector indicating the upper bounds of each parameter in \eqn{\theta}.}
#' }
#' 
#' Additional control variables can be supplied to \code{MagiSolver} via the optional list \code{control}, which may include the following:
#' \describe{
#'   \item{\code{sigma}}{a vector of noise levels (observation noise standard deviations) \eqn{\sigma} for each component, at which to initialize MCMC sampling.  By default, \code{MagiSolver} computes starting values for \code{sigma} via Gaussian process (GP) smoothing. If the noise levels are known, specify \code{sigma} together with \code{useFixedSigma = TRUE}.}
#'   \item{\code{phi}}{a matrix of GP hyper-parameters for each component, with rows for the kernel hyper-parameters and columns for the system components. By default, \code{MagiSolver} estimates \code{phi} via an optimization routine.}
#'   \item{\code{theta}}{a vector of starting values for the parameters \eqn{\theta}, at which to initialize MCMC sampling. By default, \code{MagiSolver} uses an optimization routine to obtain starting values.}
#'   \item{\code{xInit}}{a matrix of values for the system trajectories of the same dimension as \code{y}, at which to initialize MCMC sampling. Default is linear interpolation between the observed (non-missing) values of \code{y} and an optimization routine for entirely unobserved components of \code{y}.}
#'   \item{\code{mu}}{a matrix of values for the mean function of the GP prior, of the same dimension as \code{y}. Default is a zero mean function.}
#'   \item{\code{dotmu}}{a matrix of values for the derivatives of the GP prior mean function, of the same dimension as \code{y}. Default is zero.}
#'   \item{\code{priorTemperature}}{the tempering factor by which to divide the contribution of the GP prior, to control the influence of the GP prior relative to the likelihood. Default is the total number of observations divided by the total number of discretization points.}
#'   \item{\code{niterHmc}}{MCMC sampling from the posterior is carried out via the Hamiltonian Monte Carlo (HMC) algorithm. \code{niterHmc} specifies the number of HMC iterations to run.  Default is 20000 HMC iterations.}
#'   \item{\code{nstepsHmc}}{the number of leapfrog steps per HMC iteration. Default is 200.}
#'   \item{\code{burninRatio}}{the proportion of HMC iterations to be discarded as burn-in. Default is 0.5, which discards the first half of the MCMC samples.}
#'   \item{\code{stepSizeFactor}}{initial leapfrog step size factor for HMC. Can be a specified as a scalar (applied to all posterior dimensions) or a vector (with length corresponding to the dimension of the posterior). Default is 0.01, and the leapfrog step size is automatically tuned during burn-in to achieve an acceptance rate between 60-90\%.}
#'   \item{\code{bandSize}}{a band matrix approximation is used to speed up matrix operations, with default band size 20. Can be increased if \code{MagiSolver} returns an error indicating numerical instability.}
#'   \item{\code{useFixedSigma}}{logical, set to \code{TRUE} if \code{sigma} is known.  If \code{useFixedSigma = TRUE}, the known values of \eqn{\sigma} must be supplied via the \code{sigma} control variable. Default is \code{FALSE}.}
#'   \item{\code{kerneltype}}{the GP covariance kernel, \code{generalMatern} is the default and recommended choice. Other available choices are \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}. See \code{\link{calCov}} for their definitions.}
#'   \item{\code{skipMissingComponentOptimization}}{logical, set to \code{TRUE} to skip automatic optimization for missing components. If \code{skipMissingComponentOptimization = TRUE}, values for \code{xInit} and \code{phi} must be supplied for all system components. Default is \code{FALSE}.}
#'   \item{\code{positiveSystem}}{logical, set to \code{TRUE} if the system cannot be negative. Default is \code{FALSE}.}
#'   \item{\code{verbose}}{logical, set to \code{TRUE} to output diagnostic and progress messages to the console. Default is \code{FALSE}.}
#'   
#' }
#' 
#' @examples
#' # Set up odeModel list for the Fitzhugh-Nagumo equations
#' fnmodel <- list(
#'   fOde = fnmodelODE,
#'   fOdeDx = fnmodelDx,
#'   fOdeDtheta = fnmodelDtheta,
#'   thetaLowerBound = c(0, 0, 0),
#'   thetaUpperBound = c(Inf, Inf, Inf)
#' )
#' 
#' # Example noisy data observed from the FN system
#' data(FNdat)
#'
#' # Set discretization for a total of 81 equally-spaced time points from 0 to 20
#' yinput <- setDiscretization(FNdat, by = 0.25)
#' 
#' # Run MagiSolver
#' # Short sampler run for demo only, more iterations needed for convergence
#' MagiSolver(yinput, fnmodel, control = list(nstepsHmc = 5, niterHmc = 101))
#' \donttest{
#' # Use 3000 HMC iterations with 100 leapfrog steps per iteration
#' FNres <- MagiSolver(yinput, fnmodel, control = list(nstepsHmc = 100, niterHmc = 3000))
#' # Summary of parameter estimates
#' summary(FNres)
#' # Plot of inferred trajectories
#' plot(FNres, comp.names = c("V", "R"), xlab = "Time", ylab = "Level")}
#' 
#' 
#' @references 
#' Wong, S. W. K., Yang, S., & Kou, S. C. (2024). `magi`: A Package for Inference of Dynamic Systems from Noisy and Sparse Data via Manifold-Constrained Gaussian Processes. *Journal of Statistical Software*, 109 (4), 1-47. \doi{10.18637/jss.v109.i04}
#'
#' Yang, S., Wong, S. W. K., & Kou, S. C. (2021). Inference of Dynamic Systems from Noisy and 
#' Sparse Data via Manifold-constrained Gaussian Processes. *Proceedings of the National Academy of Sciences*, 118 (15), e2020397118. \doi{10.1073/pnas.2020397118}
#' 
#' 
#' @export 
MagiSolver <- function(y, odeModel, tvec, control = list()) {
  if (missing(tvec) & !("time" %in% colnames(y)))
    stop("must supply a \"time\" column in y or specify a vector of time points via the argument \"tvec\"")
  
  if (missing(tvec)) {
    tvec <- y[,"time"]
    y <- y[, !colnames(y)=="time", drop=FALSE]
  }
  
  if (!is.null(control$sigma)) {
    sigmaExogenous = control$sigma
    sigmaExogenous[!is.finite(sigmaExogenous)] <- 1
  }
  else
    sigmaExogenous = numeric(0)

  if (!is.null(control$phi))
    phiExogenous = control$phi
  else
    phiExogenous = matrix(numeric(0))

  if (!is.null(control$xInit))
    xInitExogenous = control$xInit
  else
    xInitExogenous = matrix(numeric(0))

  if (!is.null(control$theta))
    thetaInitExogenous = control$theta
  else
    thetaInitExogenous = matrix(numeric(0))

  if (!is.null(control$mu))
    muExogenous = control$mu
  else
    muExogenous = matrix(numeric(0))

  if (!is.null(control$dotmu))
    dotmuExogenous = control$dotmu
  else
    dotmuExogenous = matrix(numeric(0))

  if (!is.null(control$priorTemperature)) {
    priorTemperatureLevel = control$priorTemperature
    priorTemperatureDeriv = control$priorTemperature
  } else {
    priorTemperatureLevel = 1/mean(!is.na(y))
    priorTemperatureDeriv = 1/mean(!is.na(y))
  }

  if (!is.null(control$niterHmc))
    niterHmc = control$niterHmc
  else
    niterHmc = 20000

  if (!is.null(control$burninRatio))
    burninRatio = control$burninRatio
  else
    burninRatio = 0.5

  if (!is.null(control$nstepsHmc))
    nstepsHmc = control$nstepsHmc
  else
    nstepsHmc = 200

  if (!is.null(control$stepSizeFactor))
    stepSizeFactor = control$stepSizeFactor
  else
    stepSizeFactor = 0.01

  if (!is.null(control$bandSize))
    bandSize = control$bandSize
  else
    bandSize = 20

  if (!is.null(control$useFixedSigma))
    useFixedSigma = control$useFixedSigma
  else
    useFixedSigma = FALSE

  if (!is.null(control$skipMissingComponentOptimization))
    skipMissingComponentOptimization = control$skipMissingComponentOptimization
  else
    skipMissingComponentOptimization = FALSE

  if (!is.null(control$positiveSystem))
    positiveSystem  = control$positiveSystem 
  else
    positiveSystem  = FALSE  
    
  if (!is.null(control$verbose))
    verbose = control$verbose
  else
    verbose = FALSE

  if (!is.null(control$kerneltype))
    kernel = control$kerneltype
  else
    kernel = "generalMatern"

  samplesCpp <- solveMagiRcpp(
    yFull = data.matrix(y),
    odeModel = odeModel,
    tvecFull = tvec,
    sigmaExogenous = sigmaExogenous,
    phiExogenous = phiExogenous,
    xInitExogenous = xInitExogenous,
    thetaInitExogenous = thetaInitExogenous,
    muExogenous = muExogenous,
    dotmuExogenous = dotmuExogenous,
    priorTemperatureLevel = priorTemperatureLevel,
    priorTemperatureDeriv = priorTemperatureDeriv,
    priorTemperatureObs = 1,
    kernel = kernel,
    nstepsHmc = nstepsHmc,
    burninRatioHmc = burninRatio,
    niterHmc = niterHmc,
    stepSizeFactorHmc = stepSizeFactor,
    nEpoch = 1,
    bandSize = bandSize,
    useFrequencyBasedPrior = TRUE,
    useBand = TRUE,
    useMean = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = useFixedSigma,
    skipMissingComponentOptimization = skipMissingComponentOptimization,
    positiveSystem = positiveSystem,
    verbose = verbose)

  phiUsed <- samplesCpp$phi
  samplesCpp <- samplesCpp$llikxthetasigmaSamples

  samplesCpp <- samplesCpp[,,1]

  out <- samplesCpp[-1,1,drop=FALSE]

  llikId <- 1
  xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(y)))
  thetaId <- (max(xId)+1):(max(xId)+length(odeModel$thetaLowerBound))
  sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(y))

  burnin <- as.integer(niterHmc*burninRatio)
  
  ret <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
          xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                        dim=c(niterHmc-burnin, nrow(y), ncol(y))),
          lp=samplesCpp[llikId,-(1:burnin)],
          sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]),
          phi = phiUsed,
          odeModel = odeModel,
          tvec = tvec,
          y = data.matrix(y))

  class(ret) <- "magioutput"
  
  ret
}




#' Marginal log-likelihood for Gaussian process smoothing
#'
#' Marginal log-likelihood and gradient as a function of GP hyper-parameters phi and observation noise standard deviation sigma. For use in Gaussian process smoothing where values of phi and sigma may be optimized.
#'
#' @param phisig vector containing GP hyper-parameters phi and observation noise SD sigma. See \code{\link{calCov}} for the definitions of the hyper-parameters.
#' @param yobs vector of observations
#' @param rInput distance matrix between all time points of \code{yobs}
#' @param kerneltype the covariance kernel, types \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}, \code{generalMatern} are supported.  See \code{\link{calCov}} for their definitions.
#' 
#' @return A list with elements \code{value} and \code{grad}, which are the log-likelihood value and gradient with respect to \code{phisig}, respectively.
#' 
#' @examples 
#' # Suppose phi[1] = 0.5, phi[2] = 3, sigma = 0.1
#' gpsmoothllik(c(0.5, 3, 0.1), rnorm(10), abs(outer(0:9, t(0:9), '-')[, 1, ]))
#' 
#' @export
gpsmoothllik <- function(phisig, yobs, rInput, kerneltype = "generalMatern") {
  phisigllikC(phisig, data.matrix(yobs), rInput, kerneltype)
}
 

#' Gaussian process smoothing
#'
#' Estimate hyper-parameters \code{phi} and noise standard deviation \code{sigma} for a vector of observations using Gaussian process smoothing.
#' 
#' @param yobs vector of observations
#' @param tvec vector of time points corresponding to observations
#' @param kerneltype the covariance kernel, types \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}, \code{generalMatern} are supported.  See \code{\link{calCov}} for their definitions.
#' @param sigma the noise level (if known). By default, both \code{phi} and \code{sigma} are estimated. If a value for \code{sigma} is supplied, then \code{sigma} is held fixed at the supplied value and only \code{phi} is estimated.
#' 
#' @return A list containing the elements \code{phi} and \code{sigma} with their estimated values.
#'
#' @examples
#' # Sample data and observation times
#' tvec <- seq(0, 20, by = 0.5)
#' y <- c(-1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, 
#' 0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, 
#' -0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, 
#' -1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84)
#'
#' gpsmoothing(y, tvec)
#'  
#'
#' @export
gpsmoothing <- function(yobs, tvec, kerneltype = "generalMatern", sigma = NULL) {
  
  distInput  <- abs(outer(tvec, t(tvec),'-')[,1,])
  yInput <- data.matrix(yobs - mean(yobs))
  ret <- list()
  
  if (is.null(sigma)) {
    res <- gpsmooth(yInput, distInput, kerneltype, sigmaExogenScalar = -1, TRUE)
    ret$sigma <- tail(res, 1)
    ret$phi <- res[-length(res)]
  } else {
    res <- gpsmooth(yInput, distInput, kerneltype, sigmaExogenScalar = sigma, TRUE)
    ret$sigma <- sigma
    ret$phi <- res
  }
  
  ret
}
