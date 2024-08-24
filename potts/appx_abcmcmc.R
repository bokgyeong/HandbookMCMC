rm(list=ls())
require(coda); require(MASS)
require(Rcpp); require(RcppArmadillo)
require(doParallel); require(foreach)


repID = 1
eps_threshold = 12
cycle = 100
mm = 1000


#==============================================================================-
# load data, posterior samples, and auxiliary samples ----
#==============================================================================-

load('data/sim.RData')
load(paste0('postSamp/abcmcmcE', eps_threshold, 'm', cycle, '.RData'))
load(paste0('aux/aux', repID, 'm', mm, '.RData'))


dirn = paste0('appx/abcmcmc/E', eps_threshold, 'm', cycle)
ifelse(!dir.exists(dirn), dir.create(dirn), F)
filename = paste0(dirn, '/appx', repID, '.RData')



#==============================================================================-
# estimate u and d for each posterior sample ----
#==============================================================================-
postSamp = postSamp$postSamples
niter = length(postSamp)
nn = 300000
burn = 100000
postSamp = postSamp[(burn+1):(burn+nn)]

th = unique(postSamp)
nth = length(th)

invcovth = 1/var(postSamp)


outers = unique(c(0, seq(0, nth, by = 10000), nth))

start = 1; appx = c(); timeappx = 0
# load(filename); start = which(outers == nrow(appx))

sourceCpp("src/RcppFtns.cpp")


for(j in (start+1):length(outers) ){
  
  thmatj = th[(outers[j-1]+1):outers[j]]
  
  ptm = proc.time()[3]
  dummy =  cpp_appx(thmatj, stat, invcovth, thu, aux) 
  timeappx = timeappx + proc.time()[3] - ptm
  
  appx = rbind(appx, dummy)
  save(postSamp, th, nth, stat, appx, timeappx, file = filename)
  
  print(paste0('Completed ', nrow(appx), ' samples...'))
}


