rm(list=ls())
require(foreach)

source("src/RFtns.R")


repID = 1
Nin = 3


#==============================================================================-
# ACD ----
#==============================================================================-

load(paste0('appx/dmh/m', Nin, '/appx', repID, '.RData'))

dirn = paste0('acd/dmh/m', Nin)
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/acd', repID, '.RData')



NN = 10000
nn = 200000

p = 2
nrep = sapply(1:nth, function(j) sum(th[j,1] == postSamp[,1]))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acddmh = compACD(D[1:nn,], NN)
timeacddmh = proc.time()[3] - ptm

cat('m =', Nin, ':', acddmh, '\n')

save(acddmh, timeacddmh, file = filename)

