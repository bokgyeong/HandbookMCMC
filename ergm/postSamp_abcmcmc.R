rm(list=ls())
require(ergm); require(Rcpp); require(RcppArmadillo)


eps_threshold = 1.9
cycle = 20


#==============================================================================-
# load data ----
#==============================================================================-
load("data/flor.RData")

filename = paste0('postSamp/abcmcmcE', eps_threshold, 'm', cycle, '.RData')


#==============================================================================-
# ABC-MCMC ----
#==============================================================================-
burn = 100000
outer = burn + 200000
th = matrix(hat, 1, 2)

sourceCpp("src/RcppFtns.cpp")
ptm = proc.time()[3]
postSamp = ergmABCMCMC(X, COV, th, outer, cycle, epsilon = 50, eps_threshold)
rtime = proc.time()[3] - ptm

save(burn, outer, postSamp, rtime, file = filename)



