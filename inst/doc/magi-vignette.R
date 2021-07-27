## ---- include = FALSE---------------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)

## ----setup--------------------------------------------------------------------
library(magi)

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
fnmodelODE <- function(theta,x,t) {
  V <- x[,1]
  R <- x[,2]

  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = theta[3] * (V - V^3 / 3.0 + R)
  result[,2] = -1.0/theta[3] * ( V - theta[1] + theta[2] * R)
  
  result
}

## -----------------------------------------------------------------------------
# Gradient with respect to system components X
fnmodelDx <- function(theta,x,t) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  V = x[,1]
  
  resultDx[,1,1] = theta[3] * (1 - V^2)
  resultDx[,2,1] = theta[3]
  
  resultDx[,1,2] = (-1.0 / theta[3])
  resultDx[,2,2] = ( -1.0*theta[2]/theta[3] )
  
  resultDx
}

## -----------------------------------------------------------------------------
# Gradient with respect to parameters theta 
fnmodelDtheta <- function(theta,x,t) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  V = x[,1]
  R = x[,2]
  
  resultDtheta[,3,1] = V - V^3 / 3.0 + R
  
  resultDtheta[,1,2] =  1.0 / theta[3] 
  resultDtheta[,2,2] = -R / theta[3]
  resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)
  
  resultDtheta
}

## -----------------------------------------------------------------------------
testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDtheta, 
    "FN equations", cbind(V,R), c(.5, .6, 2), tvec)

## -----------------------------------------------------------------------------
fnmodel <- list(
  fOde=fnmodelODE,
  fOdeDx=fnmodelDx,
  fOdeDtheta=fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf)
)

## -----------------------------------------------------------------------------
yobs <- data.frame(time=tvec, V=V, R=R)  

## ---- results="hide"----------------------------------------------------------
yinput <- setDiscretization(yobs, level=1)
result <- MagiSolver(yinput, fnmodel, control=list(niterHmc=2000, nstepsHmc=100))

## ---- fig.align = "center", fig.width=6, fig.asp=1----------------------------
oldpar <- par(mfrow=c(2,2), mar=c(5,2,1,1))
theta.names <- c("a", "b", "c")
for (i in 1:3) {
	plot(result$theta[,i], main=theta.names[i], type="l", ylab="")
}
plot(result$lp, main="log-post", type="l", ylab="")

## -----------------------------------------------------------------------------
theta.est <- apply(result$theta, 2,
    function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
colnames(theta.est) <- theta.names
rownames(theta.est) <- c("Mean", "2.5%", "97.5%")
signif(theta.est, 3)

## ---- fig.align = "center", fig.width=7, fig.asp=0.45-------------------------
par(mfrow=c(1,2), mar=c(4,2,1,1))
compnames <- c("V", "R")
ylim_lower <- c(-3, -2)
ylim_upper <- c(3, 2)
times <- yinput[,1]

xLB <- apply(result$xsampled, c(2,3), function(x) quantile(x, 0.025))
xMean <- apply(result$xsampled, c(2,3), mean)
xUB <- apply(result$xsampled, c(2,3), function(x) quantile(x, 0.975))

for (i in 1:2) {
  plot(times, xMean[,i], type="n", xlab="time", ylab="", ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i])

  polygon(c(times, rev(times)), c(xUB[,i], rev(xLB[,i])),
          col = "skyblue", border = NA)
  points(times, yinput[,i+1], col = "grey50")
  
  lines(times, xMean[,i], lwd=1)
}

## ---- results="hide"----------------------------------------------------------
yinput2 <- setDiscretization(yobs, level=2)
result2 <- MagiSolver(yinput, fnmodel, control=list(niterHmc=2000, nstepsHmc=100))

## -----------------------------------------------------------------------------
theta.est <- apply(result2$theta, 2,
    function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
colnames(theta.est) <- theta.names
rownames(theta.est) <- c("Mean", "2.5%", "97.5%")
signif(theta.est, 3)

## -----------------------------------------------------------------------------
hes1modelODE <- function(theta, x, t) {
     P = x[,1]
     M = x[,2]
     H = x[,3] 
     
     PMHdt = array(0, c(nrow(x), ncol(x)))
     PMHdt[,1] = -theta[1]*P*H + theta[2]*M - theta[3]*P
     PMHdt[,2] = -theta[4]*M + theta[5]/(1+P^2)
     PMHdt[,3] = -theta[1]*P*H + theta[6]/(1+P^2) - theta[7]*H
     
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
y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5),])
names(y) <- c("time", "P", "M", "H")
y$P <- y$P * exp(rnorm(nrow(y), sd=param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd=param.true$sigma[2]))

## -----------------------------------------------------------------------------
y$H <- NaN
y$P[y$time %in% seq(7.5,240,by=15)] <- NaN
y$M[y$time %in% seq(0,240,by=15)] <- NaN

## ---- fig.align = "center", fig.width=6, fig.asp=1----------------------------
matplot(x[, "time"], x[, -1], type="l", lty=1, xlab="Time (min)", ylab="Level")
matplot(y$time, y[,-1], type="p", col=1:(ncol(y)-1), pch=20, add = TRUE)
legend("topright", c("P", "M", "H"), lty=1, col=c("black", "red", "green"))

## -----------------------------------------------------------------------------
y[,2:4] <- log(y[,2:4])

## -----------------------------------------------------------------------------
hes1logmodelODE <- function (theta, x, t) {
  eP = exp(x[, 1])
  eM = exp(x[, 2])
  eH = exp(x[, 3])
  
  PMHdt <- array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * eH + theta[2] * eM/eP - theta[3]
  PMHdt[, 2] = -theta[4] + theta[5]/(1 + eP^2)/eM
  PMHdt[, 3] = -theta[1] * eP + theta[6]/(1 + eP^2)/eH - theta[7]
  PMHdt
}

## -----------------------------------------------------------------------------
hes1logmodelDx <- function (theta, x, t) {
  P = x[, 1]
  M = x[, 2]
  H = x[, 3]
  
  Dx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  
  dP = -(1 + exp(2 * P))^(-2) * exp(2 * P) * 2
  Dx[, 1, 1] = -theta[2] * exp(M - P)
  Dx[, 2, 1] = theta[2] * exp(M - P)
  Dx[, 3, 1] = -theta[1] * exp(H)
  Dx[, 1, 2] = theta[5] * exp(-M) * dP
  Dx[, 2, 2] = -theta[5] * exp(-M)/(1 + exp(2 * P))
  Dx[, 1, 3] = -theta[1] * exp(P) + theta[6] * exp(-H) * dP
  Dx[, 3, 3] = -theta[6] * exp(-H)/(1 + exp(2 * P))
  
  Dx
}

## -----------------------------------------------------------------------------
hes1logmodelDtheta <- function (theta, x, t) {
  
  P = x[, 1]
  M = x[, 2]
  H = x[, 3]
  
  Dtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  Dtheta[, 1, 1] = -exp(H)
  Dtheta[, 2, 1] = exp(M - P)
  Dtheta[, 3, 1] = -1
  Dtheta[, 4, 2] = -1
  Dtheta[, 5, 2] = exp(-M)/(1 + exp(2 * P))
  Dtheta[, 1, 3] = -exp(P)
  Dtheta[, 6, 3] = exp(-H)/(1 + exp(2 * P))
  Dtheta[, 7, 3] = -1
  
  Dtheta
}

## -----------------------------------------------------------------------------
hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = hes1logmodelDx,
  fOdeDtheta = hes1logmodelDtheta,
  thetaLowerBound = rep(0,7),
  thetaUpperBound = rep(Inf,7)
)

## ---- results='hide'----------------------------------------------------------
resultHes1 <- MagiSolver(y, hes1logmodel,
           control=list(sigma = c(0.15,0.15,NA), useFixedSigma = TRUE))

## ---- fig.align = "center", fig.width=7.2, fig.asp=0.5------------------------
par(mfrow=c(2,4), mar=c(5,2,1,1))
theta.names <- c("a", "b", "c", "d", "e", "f", "g")
for (i in 1:7) {
  plot(resultHes1$theta[,i], main=theta.names[i], type="l", ylab="")
}
plot(resultHes1$lp, main="log-posterior", type="l", ylab="")

## -----------------------------------------------------------------------------
theta.est <- apply(resultHes1$theta, 2,
    function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
colnames(theta.est) <- theta.names
rownames(theta.est) <- c("Post.Mean", "2.5%", "97.5%")
signif(theta.est, 3)

## ---- fig.align = "center", fig.width=7.2, fig.asp=0.4------------------------
ylim_lower <- c(1.5, 0.5, 0)
ylim_upper <- c(10.0, 3.5, 21)

layout(rbind(c(1,2,3), c(4,4,4)), heights = c(5,1))
compnames <- c("P", "M", "H")
compobs <- c("17 observations", "16 observations", "unobserved")
times <- y[,1]

xLB <- exp(apply(resultHes1$xsampled, c(2,3), function(x) quantile(x, 0.025)))
xMean <- exp(apply(resultHes1$xsampled, c(2,3), mean))
xUB <- exp(apply(resultHes1$xsampled, c(2,3), function(x) quantile(x, 0.975)))

for (i in 1:3) {
  plot(times, xMean[,i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(paste0(compnames[i], " (", compobs[i], ")"), cex=1)  

  polygon(c(times, rev(times)), c(xUB[,i], rev(xLB[,i])),
          col = "skyblue", border = NA)
  
  lines(x[,1], x[,1+i], col="red", lwd=2)
  lines(times, xMean[,i], col="forestgreen", lwd=2)
}

par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)

legend("center", c("truth", "inferred trajectory", "95% interval"), lty=c(1,1,0), lwd=c(2,2,0),
       col = c("red", "forestgreen", NA), fill=c(0, 0,"skyblue"), text.width=c(0, 0.4, 0.05), bty = "n",
       border=c(0, 0, "skyblue"), pch=c(NA, NA, 15), horiz=TRUE)

## -----------------------------------------------------------------------------
theta.names <- c("lambda", "rho", "delta", "N", "c")

hivtdmodelODE <- function(theta,x,tvec) {
  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = lambda - rho * TU - eta * TU * V
  result[,2] = eta * TU * V - delta * TI
  result[,3] = N * delta * TI - c * V

  result
}

hivtdmodelDx <- function(theta,x,tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))

  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  resultDx[,1,1] = -rho - eta * V
  resultDx[,2,1] = 0
  resultDx[,3,1] = -eta * TU

  resultDx[,1,2] = eta * V
  resultDx[,2,2] = -delta
  resultDx[,3,2] = eta * TU

  resultDx[,1,3] = 0
  resultDx[,2,3] = N * delta
  resultDx[,3,3] = -c

  resultDx
}

hivtdmodelDtheta <- function(theta,x,tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))

  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  resultDtheta[,1,1] = 1
  resultDtheta[,2,1] = -TU
  resultDtheta[,3,1] = 0
  resultDtheta[,4,1] = 0
  resultDtheta[,5,1] = 0

  resultDtheta[,1,2] = 0
  resultDtheta[,2,2] = 0
  resultDtheta[,3,2] = -TI
  resultDtheta[,4,2] = 0
  resultDtheta[,5,2] = 0

  resultDtheta[,1,3] = 0
  resultDtheta[,2,3] = 0
  resultDtheta[,3,3] = N * TI
  resultDtheta[,4,3] = delta * TI
  resultDtheta[,5,3] = -V

  resultDtheta
}

## -----------------------------------------------------------------------------
param.true <- list(
  theta = c(36, 0.108, 0.5, 1000, 3), # lambda, rho, delta, N, c
  x0 = c(600, 30, 1e5), # TU, TI, V
  sigma= c(sqrt(10), sqrt(10), 10)
)

## -----------------------------------------------------------------------------
times <- seq(0, 20, 0.1)

modelODE <- function(t, state, parameters) {
  list(as.vector(hivtdmodelODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = param.true$x0, times = times,
                      func = modelODE, parms = param.true$theta)
xsim <- xtrue
set.seed(12321)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=param.true$sigma[j])
}

## ---- fig.align = "center", fig.width=6, fig.asp=0.9--------------------------
matplot(xsim[,"time"], xsim[,-1], type="p", col=1:(ncol(xsim)-1),
        pch=20, log = 'y', ylab="Concentration", xlab="time")
legend("topright", c("TU", "TI", "V"), pch=20, col=c("black", "red", "green"))

## -----------------------------------------------------------------------------
hivtdmodel <- list(
  fOde=hivtdmodelODE,
  fOdeDx=hivtdmodelDx,
  fOdeDtheta=hivtdmodelDtheta,
  thetaLowerBound=c(0,0,0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf,Inf,Inf)
)

y <- setDiscretization(data.frame(xsim), 0)

## -----------------------------------------------------------------------------
testDynamicalModel(hivtdmodelODE, hivtdmodelDx, hivtdmodelDtheta,
    "HIV time-dependent system", y[,2:4], param.true$theta, y[,"time"])

## -----------------------------------------------------------------------------
phiExogenous <- matrix(0, nrow=2, ncol=ncol(y)-1)
sigmaInit <- rep(0, ncol(y)-1)
for (j in 1:(ncol(y)-1)){
  hyperparam <- gpsmoothing(y[,j+1], y[,"time"])
  phiExogenous[,j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
}

phiExogenous
sigmaInit

## ---- results="hide"----------------------------------------------------------
phiExogenous[,3] <- c(5e7, 1)
sigmaInit[3] <- 1
HIVresult <- MagiSolver(y, hivtdmodel,
      control = list(phi=phiExogenous, sigma=sigmaInit, niterHmc=10000))

## -----------------------------------------------------------------------------
theta.est <- apply(HIVresult$theta, 2,
    function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
colnames(theta.est) <- theta.names
rownames(theta.est) <- c("Mean", "2.5%", "97.5%")
signif(theta.est, 3)

## ---- fig.align = "center", fig.width=7, fig.asp=0.45-------------------------
par(mfrow=c(1,3), mar=c(4,3,1.5,1))
compnames <- c("TU", "TI", "V")
ylim_lower <- c(100, 0, 0)
ylim_upper <- c(750, 175, 1e5)

xMean <- apply(HIVresult$xsampled, c(2,3), mean)

for (i in 1:3) {
  plot(times, xMean[,i], type="n", xlab="time", ylab="", ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i])
  points(times, xsim[,i+1], col = "grey50")
  lines(times, xMean[,i], col="forestgreen", lwd=4)
  lines(times, xtrue[,i+1], col="red", lwd=1.5)
}

par(oldpar) # reset to previous pars

