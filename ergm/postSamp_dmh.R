rm(list=ls())
require(ergm); require(Rcpp); require(RcppArmadillo)


cycle = 3

#==============================================================================-
# load data ----
#==============================================================================-
load("data/flor.RData")

filename = paste0('postSamp/dmh', cycle, '.RData')


#==============================================================================-
# DMH ----
#==============================================================================-
burn = 1000
outer = burn + 200000
th = matrix(hat, 1, 2)


sourceCpp("src/RcppFtns.cpp")
ptm = proc.time()[3]
postSamp = ergmDMH(X, COV, th, outer, cycle)[-1,]
rtime = proc.time()[3] - ptm

save(burn, outer, postSamp, rtime, file = filename)


