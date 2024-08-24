rm(list=ls())
require(ergm); require(Rcpp); require(RcppArmadillo)


eps_threshold = 12
cycle = 100


load("data/sim.RData")

dirn = 'postSamp'
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/abcmcmcE', eps_threshold, 'm', cycle, '.RData')


COV = 1
sigma2 = 1

burn = 100000
outer = burn + 300000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


# initialize parameter
beta = 1


sourceCpp("src/RcppFtns.cpp")

ptm <- proc.time()[3]
postSamp = pottsABCMCMC(
  foo, ncolor, stat, COV, beta, outer, cycle, epsilon = 300, eps_threshold, 
  updateCOV, sigma2, adaptInterval, adaptFactorExponent, adapIter)
rtime <- proc.time()[3]-ptm

save(postSamp, rtime, burn, file = filename)
