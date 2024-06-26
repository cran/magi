library(testthat)
testthat::context("testing GP kernel cpp functions")

phitrue <- list(
  generalMatern = c(2.04, 1.313)
)

xtime <- seq(0,20,0.5)

kerneltype <- "generalMatern"
testthat::test_that("check C, dCdphi for generalMatern",{
  
  egcovR <- calCov(phitrue[[kerneltype]], 
                   as.matrix(dist(xtime)),
                   -sign(outer(xtime,xtime,'-')),
                   kerneltype = kerneltype)
  dimnames(egcovR$Cdoubleprime) <- NULL
  dimnames(egcovR$C) <- NULL
  
  egcovC <- generalMaternCovRcpp(phitrue[[kerneltype]], outer(xtime,xtime,'-'))
  expect_equal(egcovC$C, egcovR$C, check.attributes = FALSE)
  expect_equal(egcovC$dCdphiCube, egcovR$dCdphiCube, check.attributes = FALSE)
  expect_equal(egcovC$Cprime, egcovR$Cprime, check.attributes = FALSE)
  expect_equal(egcovC$Cdoubleprime, egcovR$Cdoubleprime, check.attributes = FALSE)
  expect_equal(egcovC$Cinv, egcovR$Cinv, check.attributes = FALSE)
  expect_equal(egcovC$mphi, egcovR$mphi, check.attributes = FALSE)
  expect_equal(egcovC$Kinv, egcovR$Kinv, check.attributes = FALSE, tolerance = 1e-5)
 
  dim(egcovC$Sigma)
  dim(egcovC$dSigmadphiCube)
  
  #if(interactive()){
  #  mb <- microbenchmark::microbenchmark(
  #    egcovR <- calCovGeneralMatern(phitrue[[kerneltype]],
  #                                  as.matrix(dist(xtime)),
  #                                  -sign(outer(xtime,xtime,'-'))),
  #    egcovC <- generalMaternCovRcpp(phitrue[[kerneltype]], outer(xtime,xtime,'-'))
  #  )
  #  print(mb)
  #}
  
})


testthat::test_that("generalMatern analytical / numerical gradient",{
  xtime <- seq(0,2,0.1)
  testpintPhi <- rexp(2)
  delta <- 1e-4
  egcov0 <- calCov(testpintPhi, 
                   as.matrix(dist(xtime)),
                   -sign(outer(xtime,xtime,'-')),
                   kerneltype = "generalMatern")
  
  egcov0cpp <- generalMaternCovRcpp(testpintPhi, outer(xtime,xtime,'-'))
  dimnames(egcov0$C) <- NULL
  dimnames(egcov0$dCdphiCube) <- NULL
  expect_equal(egcov0cpp$C, egcov0$C, check.attributes = FALSE, tolerance=1e-6)
  expect_equal(egcov0cpp$dCdphiCube, egcov0$dCdphiCube, check.attributes = FALSE, tolerance=1e-6)
  expect_equal(egcov0cpp$Cprime, egcov0$Cprime, check.attributes = FALSE)
  expect_equal(egcov0cpp$Cdoubleprime, egcov0$Cdoubleprime, check.attributes = FALSE)
  
  egcov0cpp$dCprimedphiCube
  egcov0cpp$dCdoubleprimedphiCube
  
  
  for(it in 1:length(testpintPhi)){
    testpintPhi1 <- testpintPhi
    testpintPhi1[it] <- testpintPhi1[it] + delta
    egcov1cpp <- generalMaternCovRcpp(testpintPhi1, outer(xtime,xtime,'-'))
    expect_equal((egcov1cpp$C - egcov0cpp$C)/delta, egcov0cpp$dCdphiCube[,,it],
                 tolerance = 1e-3, scale = max(abs(egcov0cpp$dCdphiCube[,,it])))
    
    expect_equal((egcov1cpp$Cprime - egcov0cpp$Cprime)/delta, egcov0cpp$dCprimedphiCube[,,it],
                 tolerance = 1e-3, scale = max(abs(egcov0cpp$dCprimedphiCube[,,it])))
    
    expect_equal((egcov1cpp$Cdoubleprime - egcov0cpp$Cdoubleprime)/delta, egcov0cpp$dCdoubleprimedphiCube[,,it],
                 tolerance = 1e-3, scale = max(abs(egcov0cpp$dCdoubleprimedphiCube[,,it])))
    
  }
})

testthat::test_that("check Cprime and Cdoubleprime", {
  xtime <- seq(0,10e-3,1e-3)
  delta <- mean(diff(xtime))
  testpintPhi <- rexp(2)
  egcov2cpp <- generalMaternCovRcpp(testpintPhi, outer(xtime,xtime,'-'))
  egcov2 <- calCov(testpintPhi, 
                   as.matrix(dist(xtime)),
                   -sign(outer(xtime,xtime,'-')),
                   kerneltype = kerneltype,
                   noiseInjection = 0)
  dimnames(egcov2$Cprime) <- NULL
  dimnames(egcov2$Cdoubleprime) <- NULL
  dimnames(egcov2$C) <- NULL
  testthat::expect_equal( egcov2cpp$Cprime, egcov2$Cprime)
  testthat::expect_equal( egcov2cpp$Cdoubleprime, egcov2$Cdoubleprime)
  
  testthat::expect_equal( ((egcov2$Cprime[-1,1] + egcov2$Cprime[-nrow(egcov2$Cprime),1])/2), 
                          (diff(egcov2$C[,1])/delta), 
                          tolerance = 0.0001, scale=max(abs(egcov2$Cprime[-1,1])))
  testthat::expect_equal( (egcov2$Cdoubleprime[-1,1] + egcov2$Cdoubleprime[-nrow(egcov2$Cdoubleprime),1])/2, 
                          -diff(egcov2$Cprime[,1])/delta, 
                          tolerance = 0.001, scale=max(abs(egcov2$Cdoubleprime[-1,1])))
})
