rm(list=ls())
require(potts); require(Rcpp); require(RcppArmadillo)
require(DiceKriging); require(DiceDesign)
require(fields); require(foreach)


LL = 500
dd = 20


dirn = 'postSamp'
ifelse(!dir.exists(dirn), dir.create(dirn), FALSE)
filename = paste0(dirn, '/LikEmL', LL, 'd', dd, '.RData')


#==============================================================================-
# pre-computation ----
#==============================================================================-
cycle = 100
load(paste0("postSamp/dmh", cycle, ".RData"))

hat <- mean(postSamples[2e3:3e3])
th <- as.matrix( (max(postSamples[2e3:3e3])-min(postSamples[2e3:3e3]))*runif(1000) + min(postSamples[2e3:3e3]) )

save(hat, th, stat, cycle, file = 'postSamp/preEmul.RData')



## Generate auxiliary variables for each particle ----
load("data/sim.RData")
load('postSamp/preEmul.RData')
sourceCpp("src/RcppFtns.cpp")


th_ = th
ptm <- proc.time()[3]
Sample = pResponse(foo, ncolor, cycle, hat, LL)
IStime = proc.time()[3] - ptm

save(cycle, Sample, IStime, th_, file = paste0('postSamp/preEmulIS', LL, '.RData'))



#==============================================================================-
# LikEm ----
#==============================================================================-
load("data/sim.RData")
load('postSamp/preEmul.RData')
load(paste0('postSamp/preEmulIS', LL, '.RData'))
sourceCpp("src/RcppFtns.cpp")


th = as.matrix(th_[1:dd,])
y <- c()
for(i in 1:dd){
  cal = Sample%*%(th[i]-hat)
  mx = max(cal)
  y[i] = mx + log( mean( exp(cal-mx) ) )
}


burn = 1000
outer = burn + 300000
theta <- matrix(hat,1)


### unnormalized likelihood
lhX <- c()
for(i in 1:dd){ lhX[i] = stat%*%th[i] }


### full likelihood
y <- lhX - y


### calculating initial value of normalize const
ptm <- proc.time()[3]
m <- km(~ ., design = matrix(th[1:dd,]), response = matrix(y[1:dd],dd,1),covtype = "matern3_2")
GPtime = proc.time()[3] - ptm

x.point <- data.frame( t(theta[1,]) )
COV <- 0.01

### GLS bets coefficients 
beta.hat.gls <- c(coef(m)$trend1,coef(m)$trend2)
para <- rep(0, outer); para[1] <- theta

### Kriging
Sigma <- coef(m)$sd2*(1+sqrt(3)*rdist(th)/coef(m)$range)%*%exp(-sqrt(3)*rdist(th)/coef(m)$range)
InvSigma <- solve(Sigma)
Sigmacross <- coef(m)$sd2*(1+sqrt(3)*rdist(x.point,th)/coef(m)$range)*exp(-sqrt(3)*rdist(x.point,th)/coef(m)$range)
lhXZ <- c(1,as.double(x.point))%*%beta.hat.gls + Sigmacross%*%InvSigma%*%(y-cbind(1,th)%*%beta.hat.gls)

#pred.m <- predict(m,as.matrix(x.point), "UK")
#lhXZ <- ypred <- pred.m$mean

# MCMC run 
ptm <- proc.time()[3]
postSamp = pottsLikEm(outer,para,COV, lhXZ,beta.hat.gls,c(coef(m)$range,coef(m)$sd2),as.vector(th),y,stat)
EMtime <- proc.time()[3]-ptm

rtime = IStime + GPtime + EMtime

rtime
# ts.plot(postSamp)

save(postSamp, rtime, file = filename)






