rm(list=ls())
require(foreach)

source("src/RFtns.R")


repID = 1
dd = 200
cycle = 10


#==============================================================================-
# ACD ----
#==============================================================================-

load(paste0('appx/alr/d', dd, 'm', cycle, '/appx', repID, '.RData'))

dirn = paste0('acd/alr/d', dd, 'm', cycle)
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/acd', repID, '.RData')


NN = 10000
nn = 200000

p = 2
nrep = sapply(1:nth, function(j) sum(th[j,1] == postSamp[,1]))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acdalr = compACD(D[1:nn,], NN)
timeacdalr = proc.time()[3] - ptm

save(acdalr, timeacdalr, file = filename)

cat('d =', dd, ':', acdalr, '\n')

