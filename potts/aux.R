rm(list=ls())
require(Rcpp); require(RcppArmadillo)
require(doParallel); require(foreach)

source("src/RFtns.R")


repID = 1
mm = 1000
NN = 20000

filename = paste0('aux/aux', repID, 'm', mm, '.RData')

#==============================================================================-
# select particles ----
#==============================================================================-
cycle = 90
load('data/sim.RData')
load(paste0('postSamp/dmh', cycle, '.RData'))
burn = 1000

thu = as.matrix( ( max(postSamples[(burn+1):(burn+mm)]) - min(postSamples[(burn+1):(burn+mm)]) ) * runif(mm) + 
                   min(postSamples[(burn+1):(burn+mm)]) )



#==============================================================================-
# generate auxiliary variables for each particle
#==============================================================================-
burnN = 1000

sourceCpp("src/RcppFtns.cpp")


ptm = proc.time()[3]
aux = foreach(i = 1:mm) %do% {
  potts_stat2(foo, ncolor, thu[i], NN+burnN)[-(1:burnN)]
}
timeaux = proc.time()[3] - ptm

save(burn, mm, thu, NN, burnN, stat, aux, timeaux, file = filename)

