rm(list=ls())
require(foreach)
source('src/RFtns.R')



repID = 1
eps_threshold = 12
cycle = 100




#==============================================================================-
# ACD ----
#==============================================================================-
load(paste0('appx/abcmcmc/E', eps_threshold, 'm', cycle, '/appx', repID, '.RData'))

dirn = paste0('acd/abcmcmc/E', eps_threshold, 'm', cycle)
ifelse(!dir.exists(dirn), dir.create(dirn), F)
filename = paste0(dirn, '/acd', repID, '.RData')




NN = 20000
nn = 300000

p = 1
nrep = sapply(1:nth, function(j) sum(th[j] == postSamp))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acdabcmcmc = ACD(D[1:nn,], NN)
timeacdabcmcmc = proc.time()[3] - ptm

save(acdabcmcmc, timeacdabcmcmc, file = filename)

