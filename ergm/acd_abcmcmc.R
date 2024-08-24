rm(list=ls())
require(foreach)

source("src/RFtns.R")


repID = 1
eps_threshold = 1.9
cycle = 20


#==============================================================================-
# ACD ----
#==============================================================================-

load(paste0('appx/abcmcmc/E', eps_threshold, 'm', cycle, '/appx', repID, '.RData'))

path = paste0('acd/abcmcmc/E', eps_threshold, 'm', cycle)
ifelse(!dir.exists(path), dir.create(path), FALSE)
filename = paste0(path, '/acd', repID, '.RData')


NN = 10000
nn = 200000

p = 2
nrep = sapply(1:nth, function(j) sum(th[j,1] == postSamp[,1]))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acdabcmcmc = compACD(D[1:nn,], NN)
timeacdabcmcmc = proc.time()[3] - ptm

save(acdabcmcmc, timeacdabcmcmc, file = filename)


