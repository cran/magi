% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamicalSystemModels.R
\name{testDynamicalModel}
\alias{testDynamicalModel}
\title{Test dynamic system model specification}
\usage{
testDynamicalModel(modelODE, modelDx, modelDtheta, modelName, x, theta, tvec)
}
\arguments{
\item{modelODE}{function that computes the ODEs, specified with the form \eqn{f(theta, x, t)}. See examples.}

\item{modelDx}{function that computes the gradients of the ODEs with respect to the system components. See examples.}

\item{modelDtheta}{function that computes the gradients of the ODEs with respect to the parameters \eqn{\theta}. See examples.}

\item{modelName}{string giving a name for the model}

\item{x}{data matrix of system values, one column for each component, at which to test the gradients}

\item{theta}{vector of parameter values for \eqn{\theta}, at which to test the gradients}

\item{tvec}{vector of time points corresponding to the rows of \code{x}}
}
\value{
A list with elements \code{testDx} and \code{testDtheta}, each with value \code{TRUE} if the corresponding gradient check passed and \code{FALSE} if not.
}
\description{
Given functions for the ODE and its gradients (with respect to the system components and parameters), verify the correctness of the gradients using numerical differentiation.
}
\details{
Calls \code{\link[testthat]{test_that}} to test equality of the analytic and numeric gradients.
}
\examples{
# ODE system for Fitzhugh-Nagumo equations
fnmodelODE <- function(theta,x,t) {
  V <- x[,1]
  R <- x[,2]

  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = theta[3] * (V - V^3 / 3.0 + R)
  result[,2] = -1.0/theta[3] * ( V - theta[1] + theta[2] * R)
  
  result
}

# Gradient with respect to system components
fnmodelDx <- function(theta,x,t) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  V = x[,1]
  
  resultDx[,1,1] = theta[3] * (1 - V^2)
  resultDx[,2,1] = theta[3]
  
  resultDx[,1,2] = (-1.0 / theta[3])
  resultDx[,2,2] = ( -1.0*theta[2]/theta[3] )
  
  resultDx
}

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

# Example incorrect gradient with respect to parameters theta
fnmodelDthetaWrong <- function(theta,x,t) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  V = x[,1]
  R = x[,2]
  
  resultDtheta[,3,1] = V - V^3 / 3.0 - R
  
  resultDtheta[,1,2] =  1.0 / theta[3] 
  resultDtheta[,2,2] = -R / theta[3]
  resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)
  
  resultDtheta
}

# Sample data for testing gradient correctness
tvec <- seq(0, 20, by = 0.5)
V <- c(-1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, 
0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, 
-0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, 
-1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84)
R <- c(0.94, 1.22, 0.89, 0.13, 0.4, 0.04, -0.21, -0.65, -0.31, -0.65, 
 -0.72, -1.26, -0.56, -0.44, -0.63, 0.21, 1.07, 0.57, 0.85, 1.04, 
 0.92, 0.47, 0.27, 0.16, -0.41, -0.6, -0.58, -0.54, -0.59, -1.15, 
 -1.23, -0.37, -0.06, 0.16, 0.43, 0.73, 0.7, 1.37, 1.1, 0.85, 0.23)
 
# Correct gradients
testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDtheta, 
    "FN equations", cbind(V,R), c(.5, .6, 2), tvec)
    
# Incorrect theta gradient (test fails)
testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDthetaWrong, 
    "FN equations", cbind(V,R), c(.5, .6, 2), tvec)
    

}
