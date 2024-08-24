rm(list=ls())
# require(devtools)
# install_version("ergm", version = "3.8.0", repos = "http://cran.us.r-project.org")
require(ergm)
require(Bergm) # works with ergm version 3.7.1 and 3.8.0
require(Matrix)
require(rlist)
require(MASS)
require(spatstat)
require(statmod)


# ======================================================= 
# ALGORITHMS
# =======================================================
# Notations
# p: length of theta
# mu0: prior mean of theta (px1)
# Sigma0: prior covariance of theta (pxp)
# mu: mean of variational approximation
# Sigma: covariance of variational approximation
# APL_out: output containing parameters relating to adjusted Pseudolikelihood

b0 <- function(x) {log(1 + exp(x))}
b1 <- function(x) {1/(1 + exp(-x))}                # 1st derivative
b2 <- function(x) {1/(2 + exp(-x) + exp(x))}       # 2nd derivative
phi <- function(x) {exp(-0.5*x^2) / sqrt(2*pi)}    # std normal


# extract parameters from APL_out
extr_APL <- function(APL_out){
  Theta_PL <- APL_out$Theta_PL       # PLE
  response <- APL_out$response       # y_{ij}
  Theta_MLE <- APL_out$Theta_MLE     # MLE
  W <- APL_out$W                     # W (upper triangular matrix)
  predictor <- APL_out$predictor     # change in suff stat when y_{ij} is toggled from 0 to 1
  sy <- APL_out$sy                   # suff stat of obs network
  wgts <- APL_out$wgts               # no. of occurrences for each element in predictor
  logztheta <- APL_out$logztheta     # log(z(theta_MLE))
  logM <- APL_out$logM               # log(M)
  return(list(Theta_PL, response, Theta_MLE, W, predictor, sy, wgts, logztheta, logM))
}

# compute Log pseudolikelihood
# theta is a vector(1d)
logPL1 <- function(theta, response, predictor, wgts){
  t1 <- predictor %*%theta;                    # predictor(nxp), theta(px1)
  L <- sum(wgts*(response*t1 - b0(t1)))
  return(L)
}

# theta is an array(2d)
logPL2 <- function(theta, response, predictor, wgts){
  t1 <- predictor %*% theta;                        # predictor(nxp), theta(pxN)
  L <- colSums(wgts * (response*t1 - b0(t1)) )      # element wise *??????보기
  return(L)
}

# compute Log adjusted pseudolikelihood
logAPL <- function(theta, Theta_MLE, Theta_PL, W, logM, response, predictor, wgts){
  
  gtheta <- Theta_PL + W %*% (theta - Theta_MLE)
  L <- logM + logPL1(gtheta, response, predictor, wgts)
  return(L)
}


# compute logztheta (L = no. of temperatures, K = no. of simulations)
# est = 1 (get est of logztheta from ergm)
# est = 2 (get est of logztheta from importance sampling)
APL <- function(formula, est, burnin, K=1000, L=20){
  set.seed(1);
  n <- network.size(net); #34
  formula <- as.formula(formula);
  
  out1 <- ergm(formula, control=control.ergm(MCMLE.maxit=2000));
  sy <- out1$nw.stats;            ###target.stats
  Theta_MLE <- coef(out1);
  ll <- out1$mle.lik;
  
  out2 <- ergm(formula, estimate='MPLE');
  Theta_PL <- coef(out2);
  
  out3 <- ergmMPLE(formula);
  predictor <- out3$predictor; 
  response <- out3$response;
  wgts <- out3$weights;
  
  stats <- simulate(formula, nsim=K, coef=Theta_MLE, statsonly=TRUE, control=control.simulate(MCMC.burnin=burnin))
  
  # JULIA STARTS
  t1 <- as.vector(wgts * b2(predictor %*% Theta_PL))         # find W
  nHPL <- t(predictor) %*% (t1*predictor)        # negative Hessian at Theta_PL
  R1 <- chol(nHPL)
  R2 <- chol(cov(stats))
  W <- solve(R1, R2)
  
  if (est == 1){                                 # find log M
    logztheta <- (t(Theta_MLE) %*% sy) - ll
    logM <- ll - logPL1(Theta_PL, response, predictor, wgts)
  } else if(est ==2){
    temp <- temp <- seq(0, 1, length.out=L)
    logztheta <- 0.5*n %*% (n-1) %*% log(1+exp(Theta_MLE[1]))
    for(j in 1:L){
      theta1 <- c(Theta_MLE[1], temp[j]*Theta_MLE[2:length(Theta_MLE)])
      stats <- simulate(formula, nsim=K, coef=theta1, statsonly=TRUE,
                        control=control.simulate(MCMC.burnin=burnin))
      t1 = stats[,-1]*Theta_MLE[-1]/L
      t2 = max(t1)
      t1 = t1 - t2
      logztheta <- logztheta + t2 + log(sum(exp(t1))) - log(K)
    }
    logM <- (t(Theta_MLE) %*% sy) - logztheta - logPL1(Theta_PL, response, predictor, wgts)
  }
  APL_out <- list(Theta_MLE=Theta_MLE, response=response, predictor=predictor, sy=sy, wgts=wgts,Theta_PL=Theta_PL, W=W, logM=logM, logztheta=logztheta)
  return(APL_out)
}


# compute importance weighted lower bound (based on adjusted pseudolikelihood)
IWLB_compute <- function(APL_out, mu, Sigma, mu0, Sigma0, N=1000, In=50, maxiter=10, tol=1.0e-5) {
  start_time <- proc.time()
  
  # Extract values from APL_out
  Theta_PL <- APL_out$Theta_PL
  response <- APL_out$response
  Theta_MLE <- APL_out$Theta_MLE
  W <- APL_out$W
  predictor <- APL_out$predictor
  sy <- APL_out$sy
  wgts <- APL_out$wgts
  logztheta <- APL_out$logztheta
  logM <- APL_out$logM
  
  base <- logM - 0.5*log(det(Sigma0)) + 0.5*log(det(Sigma))
  p <- length(mu)
  C0 <- chol(Sigma0)               # Cholesky factor of Sigma0
  C <- chol(Sigma)                 # Cholesky factor of Sigma
  IWLBs <- rep(0, maxiter + 1)     # store values of IWLB
  
  dif <- 1; i <- 0; w <- rep(0, N)
  set.seed(1)
  while (dif > tol && i < maxiter) {
    i <- i + 1
    samples <- as.matrix(mvrnorm(N*In, mu, Sigma))
    mus <- matrix(rep(mu, c(rep(N*In, length(mu0)))), ncol=length(mu))
    mu0s <- matrix(rep(mu0, c(rep(N*In, length(mu0)))), ncol=length(mu0))
    w1 <- 0.5*colSums(solve(C, t(samples-mus))^2) - 0.5*colSums(solve(C0, t(samples-mu0s))^2)
    Theta_MLEs <- matrix(rep(Theta_MLE, c(rep(N*In, length(mu0)))), ncol=length(mu))
    samples <- Theta_PL + as.matrix(W) %*% t(samples - Theta_MLEs)        # compute g(theta)
    if(length(mu)==2){
      w1 <- w1 + logPL2(samples, response, predictor, wgts)
    } else{
      w1 <- w1 + logAPL(samples, Theta_MLE, Theta_PL, W, logM, response, predictor, wgts)
    } 
    if (i == 1) {wmax <- max(w1)}
    
    w1 <- matrix(exp(w1 - wmax), nrow = N)
    if (i == 1) {
      IWLBs[1] <- mean(log(w1[,1])) + wmax + base            # K=1 (ELBO)
      cat("K=", i, " IWLB=", round(IWLBs[i], 3))
    }
    
    w <- ((i-1)*In*w + rowSums(w1)) / (i*In)
    IWLBs[i+1] <- mean(log(w)) + wmax + base
    dif <- (IWLBs[i+1] - IWLBs[i]) / abs(IWLBs[i])
    if (is.nan(dif)==TRUE) {dif <-0.5}
    cat("K=", i*In, " dif=", round(dif, 4), " IWLB=", round(IWLBs[i+1], 3), "\n")
  }
  
  IWLBs <- IWLBs[1:(i+1)]
  Kvals <- c(1, seq(from = In, by = In, length.out = i)) #???거도 ??????
  IWLBdur <- proc.time() - start_time
  
  return(list(IWLBs = IWLBs, Kvals = Kvals, IWLBdur = IWLBdur))
}

# Laplace approximation
Laplace <- function(formula, APL_out, mu0, Sigma0) {
  start_time <- Sys.time()
  
  Theta_PL <- as.matrix(extr_APL(APL_out)[[1]])
  response <- as.matrix(extr_APL(APL_out)[[2]])
  Theta_MLE <- as.matrix(extr_APL(APL_out)[[3]])
  W <- as.matrix(extr_APL(APL_out)[[4]])
  predictor <- as.matrix(extr_APL(APL_out)[[5]])
  wgts <- as.matrix(extr_APL(APL_out)[[7]])
  logztheta <- as.matrix(extr_APL(APL_out)[[8]])
  logM <- as.matrix(extr_APL(APL_out)[[9]])
  
  p <- length(Theta_MLE)
  iSigma0 <- solve(Sigma0)    # inverse of Sigma0
  
  # log p(y,theta) using adjusted pseudolikelihood (omit constants indep of theta)
  # optimize minimizes instead of maximize so this is the negative of log p(y,theta)
  logpy_adj <- function(theta) {
    gtheta <- Theta_PL + as.matrix(W) %*% (as.matrix(theta) - as.matrix(Theta_MLE))
    t1 <- predictor %*% gtheta
    L <- 0.5 * t(theta - mu0) %*% (iSigma0 %*% (theta - mu0)) - sum(wgts * (response * t1 - b0(t1)))
    return(L)
  }
  # gradient of logpy_adj
  glogpy_adj <- function(theta) {
    gtheta <- Theta_PL + as.matrix(W) %*% (as.matrix(theta) - as.matrix(Theta_MLE))
    t1 <- predictor %*% gtheta
    t1 <- wgts * (response - b1(t1))
    g <- iSigma0 %*% (theta - mu0) - t(W) %*% (t(predictor) %*% t1)
    return(g)
  }
  # hessian of logpy_adj
  Hlogpy_adj <- function(theta) {
    gtheta <- Theta_PL + as.matrix(W) %*% (as.matrix(theta) - as.matrix(Theta_MLE))
    t1 <- predictor %*% gtheta
    t1 <- as.vector(wgts * b2(t1))
    H <- t(W) %*% ( t(t1 * predictor) %*% predictor ) %*% W + iSigma0
    return(H)
  }
  
  # find posterior mode (L-BFGS Algorithm)
  res <- optim(par=Theta_MLE, fn=logpy_adj, gr=glogpy_adj, method="L-BFGS-B")
  mu <- res$par
  Sigma <- solve(Hlogpy_adj(mu))
  Sigma <- (Sigma + t(Sigma)) / 2
  dur <- Sys.time() - start_time
  
  list_IWLB <- IWLB_compute(APL_out, mu, Sigma, mu0, Sigma0)
  
  IWLBs <- list_IWLB$IWLBs
  Kvals <- list_IWLB$Kvals
  IWLBdur <- list_IWLB$IWLBdur
  
  cat("IWLB =", round(tail(IWLBs, 1), 3), " mu =", round(mu, 2), " sd =", round(sqrt(diag(Sigma)), 2))
  
  out <- list("mu" = mu, "Sigma" = Sigma, "IWLBs" = IWLBs, "Kvals" = Kvals, "dur" = dur, "IWLBdur" = IWLBdur)
  return(out)
}

# compute expectation of B^{(r)}(m,v)
Bexp0 <- function(m, v, nodes, wts) {
  lm <- length(m)
  mhat <- numeric(lm)
  vhat <- numeric(lm)
  
  for (i in 1:lm) {
    lgr <- function(z) {                 # gradient of log g^{(r)}
      u <- v[i] * z + m[i]
      v[i] * b1(u) / b0(u) - z
    }
    
    # Find the mode of g^{(r)} using the `uniroot` function in R
    mhat[i] <- uniroot(lgr, interval = c(-10, 10))$root  #???거도 ??????
    
    u <- m[i] + v[i]*mhat[i]
    if (u < -30) {
      c <- 0
    } else {
      c <- (b1(u)^2) / (b0(u)^2) - b2(u) / b0(u)
    }
    vhat[i] <- (1 + v[i]^2 * c)^(-0.5)
  }
  z <- sqrt(2) * vhat %*% t(nodes) + mhat #???기?????? ???????????? 
  S <- b0(v * z + m) * phi(z)
  S <- sqrt(2) * vhat * colSums(wts * t(S))
  return(S)
}

# compute expectation of B^{(r)}(m,v)
Bexp1 <- function(m, v, nodes, wts){
  lm <- length(m)
  mhat <- numeric(lm)
  vhat <- numeric(lm)
  
  for (i in 1:lm) {
    lgr <- function(z) {            # gradient of log g^{(r)}
      u <- v[i]*z + m[i]
      v[i] / (1 + exp(u)) - z
    }
    
    mhat[i] <- uniroot(lgr, interval = c(-10, 10))$root  # mode of g^{(r)}
    
    u <- m[i] + v[i] * mhat[i]
    c <- (exp(-u) - 1) / ((1 + exp(u)) * (1 + exp(-u)))
    vhat[i] <- (1 + v[i]^2 * ((1 + exp(u))^(-2) - c))^(-0.5)
  }
  z <- sqrt(2) * vhat %*% t(nodes) + mhat #???기?????? ???????????? 
  S <- b1(v * z + m) * phi(z)
  S <- sqrt(2) * vhat * colSums(wts * t(S))
  return(S)
}

# compute expectation of B^{(r)}(m,v)
Bexp2 <- function(m, v, nodes, wts){
  lm <- length(m)
  mhat <- numeric(lm)
  vhat <- numeric(lm)
  
  for (i in 1:lm){
    lgr <- function(z){                                  # gradient of log g^{(r)}
      u = v[i]*z + m[i]
      v[i]*(exp(-u)-1)/(exp(-u)+1) - z
    }
    mhat[i] <- uniroot(lgr, interval = c(-10, 10))$root  # mode of g^{(r)}
    u = m[i] + v[i]*mhat[i];
    c = (1 - 4*exp(-u) + exp(-2*u))/(1 + exp(-u))^2
    vhat[i] <- (1 + v[i]^2 *((exp(-u)-1)^2/(exp(-u)+1)^2 - c))^(-0.5);
  }
  z <- as.matrix(sqrt(2) * vhat %*% t(nodes) + mhat) #???기?????? ???????????? 
  S <- b2(v * z + m) * phi(z)
  S <- sqrt(2) * vhat * colSums(wts * t(S))
}

# NCVMP lower bound
NCVMP_LB <- function(m, v, base, response, wgts, mu0, iSigma0, mu, Sigma, nodes, wts) {
  L <- sum(wgts * (response * m - Bexp0(m, v, nodes, wts))) + base -
    0.5*t(mu - mu0) %*% (iSigma0 %*% (mu - mu0)) +
    0.5*log(det(Sigma)) - 
    0.5*sum(iSigma0 * Sigma)
  return(L)
}

# NCVMP algorithm (based on adjusted pseudo likelihood)  #????????????
NCVMP <- function(formula, APL_out, mu0, Sigma0, tol = 1.0e-5) {
  start_time <- Sys.time()
  
  Theta_PL <- as.matrix(extr_APL(APL_out)[[1]])
  response <- as.vector(extr_APL(APL_out)[[2]])
  Theta_MLE <- as.matrix(extr_APL(APL_out)[[3]])
  W <- as.matrix(extr_APL(APL_out)[[4]])
  predictor <- as.matrix(extr_APL(APL_out)[[5]])
  wgts <- as.vector(extr_APL(APL_out)[[7]])
  logztheta <- as.matrix(extr_APL(APL_out)[[8]])
  logM <- as.matrix(extr_APL(APL_out)[[9]])
  
  alpha <- predictor %*% (Theta_PL - W %*% Theta_MLE)
  beta <- predictor %*% W
  
  p <- length(Theta_MLE)
  gauss_hermite_vals <- gauss.quad(20)
  nodes <- gauss_hermite_vals$nodes
  wts <- gauss_hermite_vals$weights * exp(nodes^2)
  
  iSigma0 <- solve(Sigma0)
  base <- 0.5 * p + logM - 0.5 * log(det(Sigma0))
  
  mu <- Theta_MLE
  Sigma <- diag(0.01, p)
  dif <- 1.0
  
  m <- as.vector(alpha + beta %*% mu)
  v <- as.vector(sqrt(rowSums((beta %*% Sigma) * beta)))
  
  
  it <- 0
  LBold <- NCVMP_LB(m, v, base, response, wgts, mu0, iSigma0, mu, Sigma, nodes, wts)
  
  while (dif > tol) {
    it <- it + 1
    Sigma <- iSigma0 + t(beta) %*% (wgts * Bexp2(m, v, nodes, wts) * beta)
    Sigma <- solve(Sigma)
    
    S1 <- (response - Bexp1(m, v, nodes, wts)) * beta * wgts
    mu <- mu + Sigma %*% (colSums(S1) - iSigma0 %*% (mu - mu0))
    
    m <- as.vector(alpha + beta %*% mu)
    v <- as.vector(sqrt(rowSums((beta %*% Sigma) * beta)))
    
    LBnew <- NCVMP_LB(m, v, base, response, wgts, mu0, iSigma0, mu, Sigma, nodes, wts)
    dif <- (LBnew - LBold) / abs(LBold)
    LBold <- LBnew
    cat(it, " LB=", round(LBnew, 3), " dif=", round(dif, 3), " mu=", round(mu, 2), "\n")
  }
  
  Sigma <- 0.5 * (Sigma + t(Sigma))
  duration <- Sys.time() - start_time
  result_IWLB <- IWLB_compute(APL_out, mu, Sigma, mu0, Sigma0)
  
  cat(" LB=", round(LBold, 3), " IWLB=", round(result_IWLB$IWLBs[length(result_IWLB$IWLBs)], 3), " mu=", round(mu, 2), " sd=", round(sqrt(diag(Sigma)), 2))
  out <- list(mu = mu, Sigma = Sigma, IWLBs = result_IWLB$IWLBs, Kvals = result_IWLB$Kvals, 
              dur = duration, LB = LBold, IWLBdur = result_IWLB$IWLBdur)
  
  return(out)
}


# find vech of a matrix A
vech <- function(A) {
  n <- nrow(A)
  idx <- lower.tri(matrix(1, n, n), diag=T)
  return(A[idx==1])
}
# return row and column indices of vech(A) where A is nxn matrix
rowcolvech <- function(n) {
  I <-NULL
  for(i in 1:n){
    cc<-i:n
    I<-c(I,cc)
  } # row indices
  J <- unlist(lapply(1:n, function(i) rep(i, n - i + 1)))      # column indices
  return(list(I = I, J = J))
}
# return indices of diagonal elements in vech(A) where A is nxn matrix
diagloc <- function(n) {
  idx <- cumsum(n:0) + 1
  return(idx[-length(idx)])
}
# recover a matrix A from vech(A)
revvech <- function(vechA) {
  n <- as.integer(0.5 * (sqrt(1 + 8 * length(vechA)) - 1))
  A <- matrix(0, n, n)
  idx <- lower.tri(matrix(1, n, n), diag=T)
  A[idx] <- vechA
  return(A)
}
# Simulate y from p(y|theta)
simulation <- function(theta, nsim, burnin) {
  stats <- simulate(as.formula(formula), nsim=nsim, coef=theta, statsonly=TRUE,
                    control=control.simulate(MCMC.burnin=burnin))
  return(stats)
}


# compute log p(y,theta) and gradient for SVI1
glogpy <- function(methd, theta, theta0, sy, C, mu0, iSigma0, stats0, stats, stats1, 
                   particles, nsim, burnin, RE, base){
  p <- length(theta)
  g <- - iSigma0 %*% (theta - mu0)
  L <- t(theta) %*% sy + base + (0.5* (t(theta - mu0) %*% g))
  t <- stats0 %*% (theta - theta0)
  L <- L - max(t) - log(sum(exp(t - max(t))))
  
  if (methd == 1) {               # generate a random sample at each iteration
    stats <- simulation(theta, nsim, burnin)
    g <- g + sy - colMeans(stats)
    ESS <- nsim
  } else if (methd == 2) {        # use importance sampling
    w <- stats %*%theta - stats1
    w <- exp(w - max(w))
    w <- w / sum(w)
    ESS <- 1 / sum(w^2)
    g <- g + sy - as.vector(t(w) %*% stats)
  } else if (methd == 3) {        # use importance sampling adaptively
    thetaa <- matrix(rep(theta0, c(rep(RE, p))), ncol=RE, byrow=T)
    idx <- which.min(colSums((solve(C, particles[,1:RE]-thetaa)^2))) # index of particle closest to theta0
    w <- stats[((idx-1)*nsim+1):(idx*nsim),] %*% theta - stats1[((idx-1)*nsim+1):(idx*nsim)]
    w = exp(w - max(w))
    w <- w / sum(w)
    ESS <- 1 / sum(w^2)
    if (ESS < nsim/3) {
      RE <- RE + 1
      particles[, RE] <- theta
      stats[((RE-1)*nsim+1):(RE*nsim), ] <- simulation(theta, nsim, burnin)
      stats1[((RE-1)*nsim+1):(RE*nsim)] <- stats[((RE-1)*nsim+1):(RE*nsim), ] %*% theta
      g <- g + sy - colMeans(stats[((RE-1)*nsim + 1):(RE*nsim), ])
      
    } else {
      g <- g + sy - as.vector(t(w)%*%stats[((idx-1)*nsim+1):(idx*nsim),])
    }
  }
  return(list(L = L, g = g, RE = RE, ESS = ESS, particles = particles, stats = stats, stats1 = stats1))
}



#SVI
SVI <- function(methd, formula, APL_out, mu0, Sigma0, burnin, mu=0, Sigma=0, nsim=1, tol=1e-5, K=1000, Interval=1000, maxit=10000, seed=1){
  
  start.time <- proc.time()
  
  formula <- as.formula(formula)
  
  theta0 <- as.matrix(APL_out$Theta_MLE)                 # MLE
  logz0 <- APL_out$logztheta                  # log z(theta) evaluated at MLE
  sy <- APL_out$sy                            # suff stats of obs network
  p <- length(theta0)                         # length of theta
  nz <- as.integer(0.5 * p * (p + 1))         # no. of elements in Cholesky factor
  I <- rowcolvech(p)$I;                       # row & column indices of elements in vech
  J <- rowcolvech(p)$J; 
  Cdiag <- diagloc(p)                         # location of diagonal entries
  
  set.seed(seed)
  
  stats0 <- simulation(theta0, K, burnin)        # for estimating lower bound
  
  if (length(mu) == 1) {mu <- theta0}       # initialize mu if not provided
  if (length(Sigma) == 1) {
    C <- lower.tri(0.1 * diag(p))   # initialize C if not provided
    Sigma <- C %*% t(C)
  } else {
    C <- t(Sigma)
  }
  Cp <- vech(C);
  Cp[Cdiag] <- log(diag(C));
  
  # Methods
  stats <- matrix(0, 1, 2)
  stats1 <- numeric(2)
  particles <- matrix(0, 1, 2)
  
  if(methd == 2){
    stats <- simulation(theta0, nsim, burnin)
    stats1 <- as.vector(stats %*% theta0)
  }else if(methd == 3){
    stats <- matrix(0, maxit*nsim, p)
    stats1 <- numeric(maxit*nsim)
    particles <- matrix(0, p, maxit)
    particles[,1] <- theta0
    stats[1:nsim,] <- simulation(theta0, nsim, burnin)
    stats1[1:nsim] <- stats[1:nsim,] %*% theta0
  }
  
  iSigma0 <- solve(Sigma0)                          # inverse of Sigma0
  base <- log(K) - 0.5 * log(det(Sigma0)) - logz0   # constant in lower bound
  ESS <- numeric(maxit)                             # store effective sample size
  thetas <- matrix(0, maxit, 2)
  mus <- matrix(0, maxit, p)
  
  # parameters for step size (Adam)
  be1 = 0.9;  be2 = 0.99;   alpha=0.001;    epsilon = 1e-8;
  mtmu = numeric(p);   vtmu = numeric(p);    # adam (1st & 2nd moment vector)
  mtCp = numeric(nz);  vtCp = numeric(nz);    # adam (1st & 2nd moment vector)
  
  LBavg = 0.0;                   # average of lower bound over Interval
  it = 0;                        # no. of iterations
  dif = 1.0;                     # relative increase in lower bound
  LBold = -Inf;
  RE=1;                          # no. of particles
  
  set.seed(seed);
  while((it < maxit) && (dif > tol)){
    it <- it + 1
    s <- as.matrix(rnorm(p))
    theta <- mu + C %*% s
    thetas[it, ] <- theta
    result <- glogpy(methd, theta, theta0, sy, C, mu0, iSigma0, stats0, stats, stats1, particles, nsim, burnin, RE, base)
    L <- result$L
    g <- result$g
    RE <- result$RE
    ESS[it] <- result$ESS
    particles <- result$particles 
    stats <- result$stats 
    stats1 <- result$stats1
    
    LBavg <- LBavg + (L + sum(log(diag(C))) + 0.5*sum(s^2))/Interval
    g <- g + solve(t(C), s) 
    
    mtmu <- be1*mtmu + (1-be1)*g
    vtmu <- be2*vtmu + (1-be2)*(g^2)
    mtmuhat <- mtmu / (1-be1^it)
    vtmuhat <- vtmu / (1-be2^it)
    mu <- mu + alpha * mtmuhat / (sqrt(vtmuhat) + epsilon)
    mus[it, ] <- mu
    
    Cpgrad <- g[I]*s[J]
    Cpgrad[Cdiag] <- Cpgrad[Cdiag] * diag(C)
    
    mtCp <- be1 * mtCp + (1-be1) * Cpgrad
    vtCp <- be2 * vtCp + (1-be2) * (Cpgrad^2)
    mtCphat <- mtCp / (1-be1^it)
    vtCphat <- vtCp / (1-be2^it)
    Cp <- Cp + alpha * mtCphat / (sqrt(vtCphat) + epsilon)
    
    C <- revvech(Cp)   # Assuming you've implemented a revvech function in R
    diag(C) <- exp(diag(C))
    
    if(it %% Interval == 0){
      if(it == Interval) {
        dif <- 1.0
      }else{
        dif <- (LBavg - LBold) / abs(LBold)
      }
      LBold <- LBavg
      Sigma <- C %*% t(C)
      cat("Iteration:", it, " LB=", round(LBavg, 3), " dif=", round(dif, 4))
      cat(round(mu, 2), " ", round(sqrt(diag(Sigma)), 2))
      LBavg <- 0.0
    }
  }
  ESS <- ESS[1:it]
  thetas <- thetas[1:it,]
  # particles <- particles[,1:RE]
  mus <- mus[1:it,]
  dur <- proc.time() - start.time
  
  
  out <- list(mu = mu, Sigma = Sigma,  dur = dur, thetas=thetas,
              ESS = ESS, mus=mus)
  return(out)
  
}




# ======================================================= 
# READ DATA
# ======================================================= 
#require(ergm.count)
#data(zach)
#net <- zach
data(florentine)
net <- flomarriage


# ======================================================= 
# RUN MODELS
# ======================================================= #
# burnin <- 30000
burnin <- 10000
# formula <- net  ~ edges + kstar(2)
formula <- net ~ edges + gwesp(decay=0.2, fixed=TRUE)
p <- 2

# Adjusted Pseudo Likelihood #
#formula, est, burnin, K=1000, L=20

start <- proc.time()
APL_out <- APL(formula, 1, burnin)   # est=1: est logztheta from ergm
dur <- proc.time() - start 
cat("APL computing time", dur, "\n")
cat("Theta_MLE=", round(APL_out$Theta_MLE, 2), 
    " Theta_PL=", round(APL_out$Theta_PL, 2), 
    " logztheta=", round(APL_out$logztheta, 2), 
    " logM=", round(APL_out$logM, 2), "\n")
# save(APL_out, file ='postSamp/APL_out_f.rdata')
# APL_out <- get(load('APL_out.rdata'))

APL_out$Theta_MLE
# colMeans(bergm_out$Theta) # -1.1037117 -0.1283107

# SVI algorithm
mu0 <- rep(0, p)
Sigma0 <- 100 * diag(p)

start <- proc.time()
bergm_out<-bergm(formula)
bergm_dur <- proc.time() - start
# save(bergm_out, file='postSamp/bergm_out_f.rdata')
# bergm_out <- get(load('bergm_out_f.rdata'))

mu <- colMeans(bergm_out$Theta)
Sigma <- cov(bergm_out$Theta)

#start <- proc.time()
#SVIb100_out <- SVI(3, formula, APL_out, mu0, Sigma0, burnin=30000, nsim=100, mu=mu, Sigma=Sigma, tol=1e-10,K=1000, Interval=1000, maxit=11000, seed=1)
#SVIb100_dur <- proc.time() - start

start <- proc.time()
SVIa5_out <- SVI(1, formula, APL_out, mu0, Sigma0, burnin=30000, nsim=5, mu=mu, Sigma=Sigma, tol=1e-10,K=1000, Interval=1000, maxit=11000, seed=1)
SVIa5_dur <- proc.time() - start

# save(SVIb100_out, file='postSamp/SVIb100_out_f.rdata')
# save(SVIa5_out, file='postSamp/SVIa5_out_f.rdata')




print(bergm_dur); print(SVIa5_dur); 
#print(SVIb100_dur)
########################################################################
########################### result ####################################
#######################################################################
#APL_out <- get(load('APL_out_f.rdata'))
#bergm_out <- get(load('bergm_out_f.rdata'))
#SVIa_out <- get(load('SVIa5_out_f.rdata'))
#SVIb_out <- get(load('SVIb100_out_f.rdata'))


postSamp = mvrnorm(100000, SVIa5_out$mu, SVIa5_out$Sigma)
save(postSamp, file = 'postSamp/vi.RData')




