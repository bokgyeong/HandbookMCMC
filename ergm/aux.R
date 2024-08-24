rm(list=ls())
require(Rcpp); require(RcppArmadillo)
require(foreach); require(qrng)

source("src/RFtns.R")


repID = 1 # index for replicate 
mm = 400 # the number of particles for SNIS
NN = 10000 # the number of auxiliary variables for each particle

filename = paste0('auxRep/aux', repID, 'm', mm, '.RData')



#==============================================================================-
# select particles ----
#==============================================================================-

load("data/flor.RData")
load('postSamp/dmh1.RData') # DMH sample used


## use quasi-random sampling ----
Domain = rbind(apply(postSamp[1:20000,], 2, min), apply(postSamp[1:20000,], 2, max)) 
p = ncol(postSamp)

gname = 'GenHal' # Generalized Halton sequence
thu = ghalton(mm, p)
thu[,1] <- (Domain[2,1] - Domain[1,1]) * thu[,1] + Domain[1,1]
thu[,2] <- (Domain[2,2] - Domain[1,2]) * thu[,2] + Domain[1,2]
# plot(thu)


#==============================================================================-
# generate auxiliary variables for each particle ----
#==============================================================================-
burnN = 1000 # burn-in


sourceCpp("src/RcppFtns.cpp")

ptm = proc.time()[3]
aux = foreach(i = 1:mm) %do% {
  Gibbs3(X, thu[i,], NN+burnN)[-(1:burnN),]
}
timeaux = proc.time()[3] - ptm

save(thu, burnN, stat, aux, timeaux, file = filename)

