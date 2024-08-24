rm(list=ls())
require(ergm); require(Rcpp); require(RcppArmadillo)

sourceCpp("src/RcppFtns.cpp")

#==============================================================================-
# Read data ----
#==============================================================================-
data(florentine)
data = flomarriage
X = data[,]
formula = data ~ edges + gwesp(0.25,fixed=TRUE)
stat = Summary(X)
m = ergm(formula, estimate="MPLE")
hat = m$coef
summary(m)
COV = solve(-m$hessian)



par(mfrow = c(1,1), mar = rep(0.6, 4))
plot(flomarriage, label.cex = 1,label = network.vertex.names(flomarriage))

save(X, data, stat, hat, COV, file = "data/flor.RData")