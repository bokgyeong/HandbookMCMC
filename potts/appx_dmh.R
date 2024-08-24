rm(list=ls())
require(coda); require(MASS)
require(Rcpp); require(RcppArmadillo)
require(doParallel); require(foreach)



repID = 1
Nin = 100
mm = 1000


#==============================================================================-
# load data, posterior samples, and auxiliary samples ----
#==============================================================================-

load('data/sim.RData')
load(paste0('postSamp/dmh', Nin, '.RData'))
load(paste0('aux/aux', repID, 'm', mm, '.RData'))


dirn = paste0('appx/dmh/m', Nin)
ifelse(!dir.exists(dirn), dir.create(dirn), F)
filename = paste0(dirn, '/appx', repID, '.RData')



#==============================================================================-
# estimate u and d for each posterior sample ----
#==============================================================================-
postSamp = postSamples
niter = length(postSamp)
nn = 300000
burn = 1000
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


