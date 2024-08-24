rm(list=ls())
require(potts); require(Rcpp); require(RcppArmadillo)

cycle = 100


load("data/sim.RData")

dirn = 'postSamp'
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/dmh', cycle, '.RData')


COV = 1
sigma2 = 1

burn = 1000
outer = burn + 300000
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1
outers = unique(c(0, seq(10000, outer, by = 10000), outer))


# initialize parameter
beta = rgamma(1, 1, 1)


start = 1; postSamples = c(); Accprob = 0; rtime = 0
# load(filename)
# start = which(outers == nrow(postSamples))

sourceCpp("src/RcppFtns.cpp")

for(i in (start+1):length(outers) ){
  outeri = outers[i]-outers[i-1]
  
  ptm = proc.time()[3]
  dummy = pottsDMH(foo, ncolor, stat, COV, beta, outeri, cycle, updateCOV, 
                   sigma2, adaptInterval, adaptFactorExponent, adapIter)
  rtime = rtime + proc.time()[3] - ptm
  
  postSamples = rbind(postSamples, dummy$postSamples)
  Accprob = ( Accprob * outers[i-1] + sum(dummy$accprob) ) / outers[i]
  
  nSamples = nrow(dummy$postSamples)
  beta = dummy$postSamples[nSamples]
  COV = dummy$COV
  adapIter = dummy$adapIter
  sigma2 = dummy$sigma2
  
  save(postSamples, Accprob, rtime, beta, COV, adapIter, sigma2, file = filename)
}
