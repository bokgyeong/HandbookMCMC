rm(list=ls())
require(foreach)

source("src/RFtns.R")


repID = 1


#==============================================================================-
# ACD ----
#==============================================================================-

load(paste0('appx/vi/appx', repID, '.RData'))

dirn = paste0('acd/vi')
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/acd', repID, '.RData')



NN = 10000
nn = nrow(postSamp)

p = 2
nrep = sapply(1:nth, function(j) sum(th[j,1] == postSamp[,1]))
D = sapply(1:(p*(p+1)/2), function(i) rep(appx[,p+i], nrep))

ptm = proc.time()[3]
acdvi = compACD(D[1:nn,], NN)
timeacdvi = proc.time()[3] - ptm

save(acdvi, timeacdvi, file = filename)

cat(acdvi, '\n')

