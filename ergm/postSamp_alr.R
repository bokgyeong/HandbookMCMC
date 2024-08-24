rm(list=ls())
require(Rcpp); require(RcppArmadillo);  require(fields)


dd = 200 # the number of particles
cycle = 10 # the number of Gibbs updates for auxiliary variable simulation


#==============================================================================-
# load data ----
#==============================================================================-
load("data/flor.RData")

filename = paste0('postSamp/alr', dd, 'm', cycle, '.RData')


#==============================================================================-
# ALR algorithm ----
#==============================================================================-

aux.par = matrix(0, dd, length(stat))


## obtain auxiliary parameters using Fractional DMH ----
sourceCpp("src/RcppFtns.cpp")

ptm = proc.time()[3]
postSamp = ergmFDMH(X, COV, th = matrix(hat, 1, 2), outer = 5500, cycle)[-1,]
rtime = proc.time()[3] - ptm
save(postSamp, rtime, file = paste0('postSamp/fdmh', cycle, '.RData'))

prerun = postSamp[-(1:500),] # discard 500 as burn-in
stand = matrix(0, nrow(prerun), ncol(prerun)) # standardization   
for(i in 1:ncol(prerun)){ 
  stand[,i] = (prerun[,i] - min(prerun[,i])) / (max(prerun[,i]) - min(prerun[,i])) 
}
stand = unique(stand) # only take unique components
dmat = rdist(stand) # distance mat


## choose auxiliary parameters through min max procedure ----
ind = 1; A = 1; Ac = 2:dim(stand)[1]
aux.par[1,] = stand[ind,]

ind = which.max( dmat[,A] )
A = c(A,ind)
Ac = Ac[-which(Ac==ind)]
aux.par[2,] = stand[ind,]

for(i in 3:dd){
  dummy = max( apply( dmat[,A] , 1, min )[Ac] )
  ind = which( dmat[,A] == dummy  )
  if(ind < dim(dmat)[1]){ 
    ind = ind 
  } else{ 
    ind = ind-floor( ind/dim(dmat)[1])*dim(dmat)[1] 
    if(ind == 0){ ind = ind-(floor( ind/dim(dmat)[1]) - 1)*dim(dmat)[1]  }
  }
  A = c(A,ind)
  Ac = Ac[-which(Ac==ind)]
  aux.par[i,] = stand[ind,]
}

dist.aux.par = rdist(aux.par) # distance matrix for aux.par (for standardized version)
for(i in 1:ncol(prerun)){ 
  aux.par[,i] = (max(prerun[,i]) - min(prerun[,i])) * aux.par[,i] + min(prerun[,i])
} 


## run ALR ----
burn = 1000
outer = burn + 200000
th = matrix(hat, 1, 2)
COV = cov(aux.par)

sourceCpp("src/RcppFtns.cpp")

ptm = proc.time()[3]
postSamp = ergmAtchade(outer, cycle, COV, th, aux.par, X)[-1,]
rtime = proc.time()[3] - ptm
save(burn, postSamp, rtime, file = filename)


