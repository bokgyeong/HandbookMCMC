rm(list=ls())
require(foreach)
source('src/RFtns.R')


repID = 1
LL = 500
dd = 20



#==============================================================================-
# ACD ----
#==============================================================================-
load(paste0('appx/likem/L', LL, 'd', dd, '/appx', repID, '.RData'))


dirn = paste0('acd/likem/L', LL, 'd', dd)
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/acd', repID, '.RData')


NN = 20000
nn = 300000


p = 1
nrep = sapply(1:nth, function(j) sum(th[j] == postSamp))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acdlikem = ACD(D[1:nn,], NN)
timeacdlikem = proc.time()[3] - ptm

save(acdlikem, timeacdlikem, file = filename)



