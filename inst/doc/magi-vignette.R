## ---- include = FALSE---------------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)

## ----setup--------------------------------------------------------------------
library("magi")

## -----------------------------------------------------------------------------
tvec <- seq(0, 20, by = 0.5)
V <- c(-1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, 
0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, 
-0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, 
-1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84)
R <- c(0.94, 1.22, 0.89, 0.13, 0.4, 0.04, -0.21, -0.65, -0.31, -0.65, 
 -0.72, -1.26, -0.56, -0.44, -0.63, 0.21, 1.07, 0.57, 0.85, 1.04, 
 0.92, 0.47, 0.27, 0.16, -0.41, -0.6, -0.58, -0.54, -0.59, -1.15, 
 -1.23, -0.37, -0.06, 0.16, 0.43, 0.73, 0.7, 1.37, 1.1, 0.85, 0.23)

## -----------------------------------------------------------------------------
fnmodelODE <- function(theta, x, tvec) {
  V <- x[, 1]
  R <- x[, 2]
  
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] = theta[3] * (V - V^3 / 3.0 + R)
  result[, 2] = -1.0/theta[3] * (V - theta[1] + theta[2] * R)
  
  result
}

## -----------------------------------------------------------------------------
# Gradient with respect to system components X
fnmodelDx <- function(theta, x, tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  V = x[, 1]
  
  resultDx[, 1, 1] = theta[3] * (1 - V^2)
  resultDx[, 2, 1] = theta[3]
  
  resultDx[, 1, 2] = -1.0 / theta[3]
  resultDx[, 2, 2] = -theta[2] / theta[3]
  
  resultDx
}

## -----------------------------------------------------------------------------
# Gradient with respect to parameters theta 
fnmodelDtheta <- function(theta, x, tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  V = x[, 1]
  R = x[, 2]
  
  resultDtheta[, 3, 1] = V - V^3 / 3.0 + R
  
  resultDtheta[, 1, 2] =  1.0 / theta[3] 
  resultDtheta[, 2, 2] = -R / theta[3]
  resultDtheta[, 3, 2] = 1.0 / (theta[3]^2) * (V - theta[1] + theta[2] * R)
  
  resultDtheta
}

## -----------------------------------------------------------------------------
testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDtheta, 
    "FN equations", cbind(V,R), c(.5, .6, 2), tvec)

## -----------------------------------------------------------------------------
fnmodel <- list(
  fOde = fnmodelODE,
  fOdeDx = fnmodelDx,
  fOdeDtheta = fnmodelDtheta,
  thetaLowerBound = c(0, 0, 0),
  thetaUpperBound = c(Inf, Inf, Inf)
)

## -----------------------------------------------------------------------------
yobs <- data.frame(time = tvec, V = V, R = R)  

## ---- results="hide"----------------------------------------------------------
yinput <- setDiscretization(yobs, level = 1)
result <- MagiSolver(yinput, fnmodel, control = list(niterHmc = 2000, nstepsHmc = 100))

## ---- fig.align = "center", fig.width=6, fig.asp=1----------------------------
oldpar <- par(mfrow = c(2, 2), mar = c(5, 2, 1, 1))
theta.names <- c("a", "b", "c")
for (i in 1:3) {
	plot(result$theta[, i], main = theta.names[i], type = "l", ylab = "")
}
plot(result$lp, main = "log-post", type = "l", ylab="")

## -----------------------------------------------------------------------------
summary(result, par.names = theta.names)

## ---- fig.align = "center", fig.width=7, fig.asp=0.45-------------------------
plot(result, comp.names = c("V", "R"), xlab = "Time", ylab = "Level")

## ---- results="hide"----------------------------------------------------------
yinput2 <- setDiscretization(yobs, level = 2)
result2 <- MagiSolver(yinput, fnmodel, control = list(niterHmc = 2000, nstepsHmc = 100))

## -----------------------------------------------------------------------------
summary(result2, par.names = theta.names)

## -----------------------------------------------------------------------------
hes1modelODE <- function(theta, x, tvec) {
  P = x[, 1]
  M = x[, 2]
  H = x[, 3] 
  
  PMHdt = array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * P * H + theta[2] * M - theta[3] * P
  PMHdt[, 2] = -theta[4] * M + theta[5] / (1 + P^2)
  PMHdt[, 3] = -theta[1] * P * H + theta[6] / (1 + P^2) - theta[7] * H
  
  PMHdt
}

## -----------------------------------------------------------------------------
param.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.439, 2.037, 17.904),
  sigma = c(0.15, 0.15, NA)
)

## -----------------------------------------------------------------------------
modelODE <- function(t, state, parameters) {
  list(as.vector(hes1modelODE(parameters, t(state), t)))
}

## -----------------------------------------------------------------------------
x <- deSolve::ode(y = param.true$x0, times = seq(0, 60*4, by = 0.01),
                  func = modelODE, parms = param.true$theta)

## -----------------------------------------------------------------------------
set.seed(12321)
y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5), ])
names(y) <- c("time", "P", "M", "H")
y$P <- y$P * exp(rnorm(nrow(y), sd = param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd = param.true$sigma[2]))

## -----------------------------------------------------------------------------
y$H <- NaN
y$P[y$time %in% seq(7.5, 240, by = 15)] <- NaN
y$M[y$time %in% seq(0, 240, by = 15)] <- NaN

## ---- fig.align = "center", fig.width=6, fig.asp=1----------------------------
compnames <- c("P", "M", "H")
matplot(x[, "time"], x[, -1], type = "l", lty = 1, 
        xlab = "Time (min)", ylab = "Level")
matplot(y$time, y[,-1], type = "p", col = 1:(ncol(y)-1), pch = 20, add = TRUE)
legend("topright", compnames, lty = 1, col = c("black", "red", "green"))

## -----------------------------------------------------------------------------
y[, names(y) != "time"] <- log(y[, names(y) != "time"])

## -----------------------------------------------------------------------------
hes1logmodelODE <- function (theta, x, tvec) {
	P = exp(x[, 1])
	M = exp(x[, 2])
	H = exp(x[, 3])
	
	PMHdt <- array(0, c(nrow(x), ncol(x)))
	PMHdt[, 1] = -theta[1] * H + theta[2] * M / P - theta[3]
	PMHdt[, 2] = -theta[4] + theta[5] / (1 + P^2) / M
	PMHdt[, 3] = -theta[1] * P + theta[6] / (1 + P^2) / H - theta[7]

	PMHdt
}

## -----------------------------------------------------------------------------
hes1logmodelDx <- function (theta, x, tvec) {
  logP = x[, 1]
  logM = x[, 2]
  logH = x[, 3]
  
  Dx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  dP = -(1 + exp(2 * logP))^(-2) * exp(2 * logP) * 2
  Dx[, 1, 1] = -theta[2] * exp(logM - logP)
  Dx[, 2, 1] = theta[2] * exp(logM - logP)
  Dx[, 3, 1] = -theta[1] * exp(logH)
  Dx[, 1, 2] = theta[5] * exp(-logM) * dP
  Dx[, 2, 2] = -theta[5] * exp(-logM) / (1 + exp(2 * logP))
  Dx[, 1, 3] = -theta[1] * exp(logP) + theta[6] * exp(-logH) * dP
  Dx[, 3, 3] = -theta[6] * exp(-logH) / (1 + exp(2 * logP))
  	
  Dx
}

## -----------------------------------------------------------------------------
hes1logmodelDtheta <- function (theta, x, tvec) {
  logP = x[, 1]
  logM = x[, 2]
  logH = x[, 3]
  
  Dtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  Dtheta[, 1, 1] = -exp(logH)
  Dtheta[, 2, 1] = exp(logM - logP)
  Dtheta[, 3, 1] = -1
  Dtheta[, 4, 2] = -1
  Dtheta[, 5, 2] = exp(-logM) / (1 + exp(2 * logP))
  Dtheta[, 1, 3] = -exp(logP)
  Dtheta[, 6, 3] = exp(-logH) / (1 + exp(2 * logP))
  Dtheta[, 7, 3] = -1
  	
  Dtheta
}

## -----------------------------------------------------------------------------
hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = hes1logmodelDx,
  fOdeDtheta = hes1logmodelDtheta,
  thetaLowerBound = rep(0, 7),
  thetaUpperBound = rep(Inf, 7)
)

## ---- results='hide'----------------------------------------------------------
hes1result <- MagiSolver(y, hes1logmodel, 
                         control = list(sigma = param.true$sigma, useFixedSigma = TRUE))

## ---- fig.align = "center", fig.width=7.2, fig.asp=0.5------------------------
par(mfrow = c(2, 4), mar = c(5, 2, 1, 1))
theta.names <- c("a", "b", "c", "d", "e", "f", "g")
for (i in 1:7) {
	plot(hes1result$theta[, i], main = theta.names[i], type = "l", ylab="")
}
plot(hes1result$lp, main = "log-post", type = "l", ylab = "")

## -----------------------------------------------------------------------------
summary(hes1result, par.names = theta.names)

## ---- fig.align = "center", fig.width=7.2, fig.asp=0.4------------------------
xLB <- exp(apply(hes1result$xsampled, c(2,3), function(x) quantile(x, 0.025)))
xMean <- exp(apply(hes1result$xsampled, c(2,3), mean))
xUB <- exp(apply(hes1result$xsampled, c(2,3), function(x) quantile(x, 0.975)))

layout(rbind(c(1, 2, 3), c(4, 4, 4)), heights = c(5, 0.5))
ylim_lower <- c(1.5, 0.5, 0)
ylim_upper <- c(10.0, 3.5, 21)
compobs <- c("17 observations", "16 observations", "unobserved")
times <- y[, "time"]

for (i in 1:3) {
  plot(times, xMean[, i], type = "n", xlab = "time", ylab = compnames[i],
        ylim = c(ylim_lower[i], ylim_upper[i]))
  mtext(paste0(compnames[i], " (", compobs[i], ")"), cex = 1)  
  
  polygon(c(times, rev(times)), c(xUB[, i], rev(xLB[, i])),
          col = "skyblue", border = NA)
  
  lines(x[, 1], x[, 1+i], col = "red", lwd = 2)
  lines(times, xMean[, i], col = "forestgreen", lwd = 2)
  points(times, exp(y[, 1+i]))
}

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n', xlab = NA, ylab = NA, frame.plot = FALSE)

legend("center", c("truth", "inferred trajectory", "95% credible interval", "noisy observations"),
       lty = c(1, 1, 0, 0), lwd = c(2, 2, 0, 1), col = c("red", "forestgreen", NA, "black"),
       fill = c(0, 0, "skyblue", 0), text.width = c(0.02, 0.25, 0.05, 0.15), bty = "n",
       border = c(0, 0, "skyblue", 0), pch = c(NA, NA, 15, 1), horiz = TRUE)

## -----------------------------------------------------------------------------
hivtdmodelODE <- function(theta, x, tvec) {
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]
  
  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))
  
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] = lambda - rho * TU - eta * TU * V
  result[, 2] = eta * TU * V - delta * TI
  result[, 3] = N * delta * TI - c * V
  
  result
}

hivtdmodelDx <- function(theta, x, tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]
  
  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))
  
  resultDx[, , 1] = cbind(-rho - eta * V, 0, -eta * TU)
  resultDx[, , 2] = cbind(eta * V, -delta, eta * TU)
  resultDx[, , 3] = cbind(rep(0, nrow(x)), N * delta, -c)

  resultDx
}

hivtdmodelDtheta <- function(theta, x, tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  delta <- theta[3]
  N <- theta[4]

  resultDtheta[, , 1] = cbind(1, -TU, 0, 0, 0)
  resultDtheta[, , 2] = cbind(0, 0, -TI, 0, 0)
  resultDtheta[, , 3] = cbind(0, 0, N * TI, delta * TI, -V)

  resultDtheta
}

## -----------------------------------------------------------------------------
param.true <- list(
  theta = c(36, 0.108, 0.5, 1000, 3),
  x0 = c(600, 30, 1e5), 
  sigma = c(sqrt(10), sqrt(10), 10), 
  times = seq(0, 20, 0.2) 
)

## -----------------------------------------------------------------------------
set.seed(12321)
modelODE <- function(tvec, state, parameters) {
  list(as.vector(hivtdmodelODE(parameters, t(state), tvec)))
}

xtrue <- deSolve::ode(y = param.true$x0, times = param.true$times,
                      func = modelODE, parms = param.true$theta)
y <- data.frame(xtrue)
for(j in 1:(ncol(y) - 1)){
  y[, 1+j] <- y[, 1+j] + rnorm(nrow(y), sd = param.true$sigma[j])
}

## ---- fig.align = "center", fig.width=6, fig.asp=0.45-------------------------
compnames <- c("TU", "TI", "V")
complabels <- c("Concentration", "Concentration", "Load")
par(mfrow = c(1, 3), mar = c(4, 4, 1.5, 1))
for (i in 1:3) {
  plot(param.true$times, y[, i+1], xlab = "Time", ylab = complabels[i])
  mtext(compnames[i])
  lines(xtrue[, "time"], xtrue[, i+1], col = "red", lwd = 2)
}

## -----------------------------------------------------------------------------
hivtdmodel <- list(
  fOde = hivtdmodelODE,
  fOdeDx = hivtdmodelDx,
  fOdeDtheta = hivtdmodelDtheta,
  thetaLowerBound = rep(0, 5),
  thetaUpperBound = rep(Inf, 5)
)

y_I <- setDiscretization(y, level = 1)

## -----------------------------------------------------------------------------
testDynamicalModel(hivtdmodelODE, hivtdmodelDx, hivtdmodelDtheta,
    "HIV time-dependent system", y[, 2:4], param.true$theta, y[, "time"])

## -----------------------------------------------------------------------------
phiEst <- matrix(0, nrow = 2, ncol = ncol(y) - 1)
sigmaInit <- rep(0, ncol(y) - 1)
for (j in 1:(ncol(y) - 1)){
  hyperparam <- gpsmoothing(y[, j+1], y[, "time"])
  phiEst[, j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
}

colnames(phiEst) <- compnames

phiEst
sigmaInit

## ---- results="hide"----------------------------------------------------------
phiEst[, 3] <- c(1e7, 0.5)
sigmaInit[3] <- 100
HIVresult <- MagiSolver(y_I, hivtdmodel,
                  control = list(phi = phiEst, sigma = sigmaInit))

## -----------------------------------------------------------------------------
summary(HIVresult, par.names = c("lambda", "rho", "delta", "N", "c"))

## ---- fig.align = "center", fig.width=7, fig.asp=0.45-------------------------
par(mfrow = c(1, 3), mar = c(4, 3, 1.5, 1))
compnames <- c("TU", "TI", "V")
ylim_lower <- c(100, 0, 0)
ylim_upper <- c(750, 175, 1e5)

xMean <- apply(HIVresult$xsampled, c(2, 3), mean)

for (i in 1:3) {
  plot(y_I[, 1], xMean[, i], type = "n", xlab = "time", ylab = "",
       ylim = c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i])
  points(y_I[, 1], y_I[, i+1], col = "grey50")
  lines(y_I[, 1], xMean[, i], col = "forestgreen", lwd = 4)
  lines(xtrue[, "time"], xtrue[, i+1], col = "red", lwd = 1.5)
}

par(oldpar) # reset to previous pars

