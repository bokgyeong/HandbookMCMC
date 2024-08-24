// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]


#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>
// #define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;




// // [[Rcpp::export]]
// mat potts_stat(RawVector foo, vec theta, int cycle){
//   
//   // Obtaining namespace of potts package
//   Environment pkg = Environment::namespace_env("potts");
//   
//   // Picking up calc_t_full() function from potts package
//   Function f = pkg["potts"];
//   
//   List result = f(foo, theta, cycle);
//   mat stat = as<arma::mat>( result[9] );
//   return stat.row(cycle - 1);
// }
// 
// 
// 
// // [[Rcpp::export]]
// mat potts_stat2(RawVector foo, vec theta, int cycle){
//   
//   // Obtaining namespace of potts package
//   Environment pkg = Environment::namespace_env("potts");
//   
//   // Picking up calc_t_full() function from potts package
//   Function f = pkg["potts"];
//   
//   List result = f(foo, theta, cycle);
//   return as<arma::mat>( result[9] );
// }
// 
// 
// 
// // [[Rcpp::export]]
// List pottsDMH(RawVector foo, mat stat, mat COV, vec theta, int outer, int cycle, 
//               bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
// 
//   double logprob, u;
//   int p = theta.size();
//   vec thetaprev = zeros(p), thetaprop = zeros(p);
//   mat postSamples = zeros(outer, p), statprop = zeros(1, p);
//   double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
//   vec accprob = zeros(outer);
//   mat cholCOV = trans( chol( sigma2 * ( COV + 0.00001 * diagmat(ones(p)) ) ) );
//   
// 
//   for(int l = 0; l< outer; l++){
//     
//     if(updateCOV){
//       if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
//         rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
//         gamma1 = 1 / pow(adapIter, c1);
//         gamma2 = c0 * gamma1;
//         sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
//         
//         COV = COV + gamma1 * ( cov( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
//         cholCOV = trans( chol( sigma2 * ( COV + 0.00001 * diagmat(ones(p)) ) ) );
//         adapIter = adapIter + 1;
//       }
//     }
//     thetaprev = theta;
//     thetaprop = thetaprev + trans( cholCOV ) * randn(p);
//     
//     if(thetaprop[p-1] > 0){
//       statprop = potts_stat(foo, thetaprop, cycle);
//       
//       logprob = sum( -0.05 * trans(thetaprop) * thetaprop + 0.05 * trans(thetaprev) * thetaprev + 
//         (statprop - stat) * (thetaprev - thetaprop) );
//       
//       u = log( randu() );
//       if( u < logprob ){
//         theta = thetaprop;
//         accprob[l] = 1;
//       } 
//     }
//     postSamples.row(l) = trans(theta);
//     
//     if ( (l+1) % 100 == 0 ) {
//       Rprintf("Generated %d samples...\n", l+1);
//     } 
//   }
// 
//   return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
//                             Rcpp::Named("accprob") = accprob,
//                             Rcpp::Named("adapIter") = adapIter,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("COV") = COV);
// }




// [[Rcpp::export]]
double potts_stat(RawVector foo, int ncolor, double beta, int cycle){
  
  // Obtaining namespace of potts package
  Environment pkg = Environment::namespace_env("potts");
  
  // Picking up calc_t_full() function from potts package
  Function f = pkg["potts"];
  
  vec theta = zeros(ncolor+1);
  theta[ncolor] = beta;
  List result = f(foo, theta, cycle);
  NumericMatrix stat = result[9];
  return stat(cycle - 1, ncolor);
}



// [[Rcpp::export]]
vec potts_stat2(RawVector foo, int ncolor, double beta, int cycle){
  
  // Obtaining namespace of potts package
  Environment pkg = Environment::namespace_env("potts");
  
  // Picking up calc_t_full() function from potts package
  Function f = pkg["potts"];
  
  vec theta = zeros(ncolor+1);
  theta[ncolor] = beta;
  List result = f(foo, theta, cycle);
  mat stat = as<arma::mat>( result[9] );
  return stat.col(ncolor);
}




// [[Rcpp::export]]
List pottsDMH(RawVector foo, int ncolor, double stat, double COV, double beta, int outer, int cycle,
              bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
  
  double logprob, u, betaprev, betaprop, statprop;
  vec postSamples = zeros(outer), accprob = zeros(outer);
  double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double cholCOV = sqrt( sigma2 * COV );;
  
  
  for(int l = 0; l< outer; l++){
    
    if(updateCOV){
      if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
        rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
        gamma1 = 1 / pow(adapIter, c1);
        gamma2 = c0 * gamma1;
        sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
        
        COV = COV + gamma1 * ( var( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
        cholCOV = sqrt( sigma2 * ( COV + 0.00001 ) );
        adapIter = adapIter + 1;
      }
    }
    betaprev = beta;
    betaprop = betaprev + cholCOV * randn();
    
    if(betaprop > 0){
      statprop = potts_stat(foo, ncolor, betaprop, cycle);
      
      logprob = -0.05 * betaprop * betaprop + 0.05 * betaprev * betaprev +
        (statprop - stat) * (betaprev - betaprop);
      
      u = log( randu() );
      if( u < logprob ){
        beta = betaprop;
        accprob[l] = 1;
      }
    }
    postSamples[l] =  beta;
    
    if ( (l+1) % 100 == 0 ) {
      Rprintf("Generated %d samples...\n", l+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
                            Rcpp::Named("accprob") = accprob,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV);
}




// [[Rcpp::export]]
List pottsFDMH(RawVector foo, int ncolor, double stat, double COV, double beta, int outer, int cycle,
               bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
  
  double logprob, u, betaprev, betaprop, statprop;
  vec postSamples = zeros(outer), accprob = zeros(outer);
  double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double cholCOV = sqrt( sigma2 * COV );;
  
  
  for(int l = 0; l< outer; l++){
    
    if(updateCOV){
      if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
        rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
        gamma1 = 1 / pow(adapIter, c1);
        gamma2 = c0 * gamma1;
        sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
        
        COV = COV + gamma1 * ( var( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
        cholCOV = sqrt( sigma2 * ( COV + 0.00001 ) );
        adapIter = adapIter + 1;
      }
    }
    betaprev = beta;
    betaprop = betaprev + cholCOV * randn();
    
    if(betaprop > 0){
      statprop = potts_stat(foo, ncolor, betaprop, cycle);
      
      logprob = -0.05 * betaprop * betaprop + 0.05 * betaprev * betaprev +
        (statprop - stat) * (betaprev - betaprop);
      logprob = logprob * 0.5;
      
      u = log( randu() );
      if( u < logprob ){
        beta = betaprop;
        accprob[l] = 1;
      }
    }
    postSamples[l] =  beta;
    
    if ( (l+1) % 100 == 0 ) {
      Rprintf("Generated %d samples...\n", l+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
                            Rcpp::Named("accprob") = accprob,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV);
}




// [[Rcpp::export]]
List pottsABCMCMC(
    RawVector foo, int ncolor, double stat, double COV, double beta, int outer, int cycle, double epsilon, double eps_threshold,
    bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
  
  double logprob, u, betaprev, betaprop, statprop, distl;
  vec postSamples = zeros(outer), accprob = zeros(outer);
  double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double cholCOV = sqrt( sigma2 * COV );
  int percentile, count = 0, nAdap = 200;
  vec dummyvec;
  vec epsilonvec(outer);
  mat dist = epsilon * ones(1,1), dummymat = ones(1,1);
  
  
  for(int l = 0; l< outer; l++){
    
    if(updateCOV){
      if( (l+1 >= adaptInterval) && (l+1 - (adaptInterval * trunc((l+1) / adaptInterval)) == 0) ){
        rhat = sum( accprob.rows(l+1-adaptInterval, l-1) ) / (adaptInterval-1);
        gamma1 = 1 / pow(adapIter, c1);
        gamma2 = c0 * gamma1;
        sigma2 = exp( log(sigma2) + gamma2 * (rhat - ropt) );
        
        COV = COV + gamma1 * ( var( postSamples.rows(l+1-adaptInterval, l-1) ) - COV );
        cholCOV = sqrt( sigma2 * COV );
        adapIter = adapIter + 1;
      }
    }
    
    if ( (dist.size() == nAdap) && (epsilon > eps_threshold) ) {
      percentile = 0.9 * dist.size();
      dummyvec = sort(dist);
      epsilon = dummyvec[percentile];
      
      dist = epsilon * ones(1,1);
      count = 0; 
    }
    
    betaprev = beta;
    betaprop = betaprev + cholCOV * randn();
    
    if(betaprop > 0){
      statprop = potts_stat(foo, ncolor, betaprop, cycle);
      
      distl = sqrt( sum((stat - statprop) * (stat - statprop)) );
      
      if( distl < epsilon ){
        logprob = -0.05 * betaprop * betaprop + 0.05 * betaprev * betaprev;  
        
        u = log( randu() );
        if( u < logprob ){
          beta = betaprop;
          accprob[l] = 1;
          
          if(count == 0){ dist(0,0) = distl; } 
          else {
            dummymat(0,0) = distl;
            dist.insert_rows(count, dummymat); 
          }
          count = count + 1;
        }
      }
    }
    postSamples[l] =  beta;
    epsilonvec[l] = epsilon;
    
    if ( (l+1) % 100 == 0 ) {
      Rprintf("Generated %d samples...\n", l+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postSamples") = postSamples,
                            Rcpp::Named("accprob") = accprob,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV,
                            Named("count") = count,
                            Named("dist") = dist,
                            Named("epsilonvec") = epsilonvec);
}





// [[Rcpp::export]]
// Atchade's adpative algorithm for Potts model
vec pottsALR(
    RawVector foo, int ncolor, double stat,
    int outer, int inner, double initial, double sigma, vec th){
  
  int d = th.n_rows, Vmin;
  vec Vis= zeros(d), gamVec(4), cVec;    // gamVec length and sequence both can be changed
  for(int i = 0; i< 4; i++){ gamVec[i] = pow( 0.1, i ); }
  
  vec parameter(outer);
  RawVector Data = foo;
  int Stat0 = stat, Stat;
  parameter(0) = initial;
  mat Esum = zeros(outer,d); // summary statistics will be stored
  
  // approximate cVec ( logZ(theta_{i}) ) until gam is small
  for(int igam = 0; igam< 4; igam++){
    
    if( igam == 3 ){ Vmin=1000; }else{ Vmin=0; } //Vmin=1000 can be changed
    double gam = gamVec(igam);
    foo = Data, Stat = Stat0;
    vec VisTemp = zeros(d);
    
    // pos denotes the variable I in Algorithm 2.1
    int pos = d-1;
    
    // cVec is the variable c in Algorithm 2.1
    if( igam == 0 ){ cVec = zeros(d); }
    
    // Stopping criteria
    while( VisTemp.min() <= Vmin || abs(  VisTemp-mean(VisTemp)*ones(d)  ).max() >= 0.2*mean(VisTemp) ){
      // Update X_{n}
      // X = MH(X, th(pos), inner);
      Stat = potts_stat(foo, ncolor, th(pos), inner);
      
      // Update I; meaning pos
      double mx = (th*Stat - cVec).max();
      vec A = exp((th*Stat - cVec - mx)) / sum( exp((th*Stat - cVec - mx))	);
      double u = randu(), l = -1, om = 0;
      while( om<u ){
        l = l+1;
        om = om +A[l];
      }
      pos = l;
      
      // Update c
      cVec = cVec + log(1+gam)*A;
      
      // We need the next two updates to control the sequence gamma
      VisTemp[pos] = VisTemp[pos] + 1;
      //if( igam ==3 ){
      //Vis[pos] = Vis[pos] + 1;
      //Esum(Vis[pos]-1,pos) = Stat;
      //}
      
    }
  }
  
  // We are now ready the start the MCMC
  // From here, gamma is deterministic small value (0.001)
  // Initilization of ${theta_n},X,I,c$
  double theta, thetaprop;
  theta =  parameter(0);
  foo = Data, Stat = Stat0;
  int pos = d-1;
  double gam = gamVec(3),prob;
  
  // Proposal Parameters
  double b = -1.5, Acc = 0, tau = 0.3;
  
  // c = cVec from above iteration
  // Bandwidth for smoothing the estimate
  int hpar = 10;
  for(int k = 0; k< outer-1; k++){
    
    // Update X_{n}
    Stat = potts_stat(foo, ncolor, th(pos), inner);
    
    // Update I; meaning pos
    double mx = (Stat*th - cVec).max();
    vec A = exp((Stat*th - cVec - mx)) / sum( exp((Stat*th - cVec - mx))	);
    double u = randu(), l = -1, om = 0;
    while( om<u ){
      l = l+1;
      om = om +A[l];
    }
    pos = l;
    
    // Update c
    cVec = cVec + log(1+gam)*A;
    Vis[pos] = Vis[pos] + 1;
    Esum(Vis[pos]-1,pos) = Stat;
    
    // Update theta if there are enough observations
    // Evaluate \pi(theta)
    // Distance to each theta^i
    uvec ind = sort_index(  abs( theta*ones(d) - th )    );
    vec cw = (1/hpar)*ones(hpar);
    
    // There might be some trouble here if none of the closest 'hpar'
    // particles around thetaprop has received data
    vec expEtheta = zeros(hpar);
    for(int i = 0; i< hpar; i++){
      if( Vis[ind[i]] == 0 ){ expEtheta[i] = 0; }else{
        vec tmp = ( theta - th[ind[i]] )*Esum.col( ind[i] ) ;
        tmp = tmp(  span( 0, Vis[ind[i]]-1 )  );
        double tmpMax = tmp.max();
        expEtheta[i] = tmpMax + log( sum( exp( tmp - tmpMax ) ) ) - log( Vis[ind[i]] ); // eq(8)'s [   ] part
      }
    }
    vec dummy2(hpar);
    for(int i = 0; i< hpar; i++){ dummy2[i] = cVec[ ind( span(0,hpar-1) )[i] ] + expEtheta[i] + cw[i];	}
    double ftheta = dummy2.max();
    // Eth = log(exp(E(x,theta))) - log(Z(theta))
    // Stat0*theta: log(exp(E(x,theta)))= E(x,theta)
    // log(sum( exp(cVec[ind[1:hpar]]+expEtheta+cw) )): log(Z(theta))
    double Eth = Stat0*theta - ftheta - log(sum( exp(dummy2 - ftheta) ));
    
    // Propose a new theta
    thetaprop = theta + sigma*randn();
    if( thetaprop <= 0 ){ prob = 0; }else{
      
      // Evaluate posterior at the new theta
      ind = sort_index(  abs( thetaprop*ones(d) - th )    );
      cw = (1/hpar)*ones(hpar);
      
      // There might be some trouble here if none of the closest 'hpar'
      // particles around thetaprop has received data
      vec expEthetaprop = zeros(hpar);
      for(int i = 0; i< hpar; i++){
        if( Vis[ind[i]] == 0 ){ expEthetaprop[i] = 0; }else{
          vec tmp = ( thetaprop - th[ind[i]] )*Esum.col( ind[i] ) ;
          tmp = tmp(  span( 0, Vis[ind[i]]-1 )  );
          double tmpMax = tmp.max();
          expEthetaprop[i] = tmpMax + log( sum( exp( tmp - tmpMax ) ) ) - log( Vis[ind[i]] ); // eq(8)'s [   ] part
        }
      }
      for(int i = 0; i< hpar; i++){ dummy2[i] = cVec[ ind( span(0,hpar-1) )[i] ] + expEthetaprop[i] + cw[i];	}
      double fthetaprop = dummy2.max();
      double Ethprop = Stat0*thetaprop - fthetaprop - log(sum( exp(dummy2 - fthetaprop) ));
      
      // Acceptance prob.
      prob = exp(Ethprop-Eth);
    }
    // Accept - Reject
    u = randu();
    if(u <= prob){ theta = thetaprop; }
    
    // Adaptive scaling of the proposal
    vec MIN(2);
    MIN[0] = 1, MIN[1] = prob;
    Acc = Acc + (1/(k+1))*( MIN.min() - Acc );
    b = b +(1/(k+1))*( MIN.min() - tau );
    
    // Save output
    parameter(k+1) = theta;
  }
  
  
  
  return(parameter);
}




// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles
mat pAuxSamp(RawVector foo, int ncolor, int cycle, vec Designmat, int m, int num){
  
  int thnrow = Designmat.n_elem;             // number of design points   
  mat H(thnrow,m);                       // cube summary statistics will be stored. (thnrow by m)
  omp_set_num_threads(num);
  
  
  //////    START OF BIGGEST CHAIN (M)     ////// 
  for(int M = 0; M < m; M++){               // m is number of importance sampling estimate 
    
    int i;
#pragma omp parallel shared(Designmat) private(i)
{	
#pragma omp for schedule(static)  
  for(i = 0; i < thnrow; i++){
    double thetaprop = Designmat(i);        
    double sumstat = potts_stat(foo, ncolor, thetaprop, cycle);
    H(i,M) = sumstat;	 	
  }
}
  }
  return(H);        	
}



// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles to construct IS estimate
vec pResponse(RawVector foo, int ncolor,  int cycle, double hatparameter, int m){
  
  vec H(m);                         // matrix where summary statistics will be stored. (m by 2)
  for(int M = 0; M < m; M++){                            // m is number of importance sampling estimate 
    // summary statistics
    H(M) = potts_stat(foo, ncolor, hatparameter, cycle);
  }
  return(H);        	
}




// [[Rcpp::export]]
// LikEm for Potts
vec pottsLikEm(int Niter, vec theta, double COV, double lhXZ, vec betahat, vec phihat, vec Designmat, vec y, double stat){
  int thnrow = Designmat.n_elem;                                        // number of design points
  double thetaprev;                                                   // befor propose in MCMC, previous parameters
  double lhXZp,logprob,u;                                               // used in MCMC step
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], sigmasqhat = phihat[1]; 
  
  mat h1(thnrow,thnrow);   // used to construct Sigma in Gaussian process
  vec h1dcross(thnrow);     // used to construct cross Sigma in Gaussian process
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i)-Designmat(j));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat);
  mat InvSigma = inv(Sigma);	    // Inverse of Sigma in Gaussian process
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	// make design matrix for linear model 
  
  
  // Start of MCMC Chain 
  for(int k = 0; k< Niter-1; k++){
    
    double Znormal = randn(); // multivariate proposal by using Cholesky factorization
    thetaprev = theta(k);
    
    // proposed parameter and corresponding rho coefficients
    double thetaprop = thetaprev + Znormal*COV;
    
    // constranits on prior space
    if( thetaprop < 0 ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  // Caculating cross covaraince matrix
        h1dcross[i] =  fabs(thetaprop-Designmat(i));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat);
      vec xpoint = ones(1);
      xpoint.insert_rows(1,thetaprop);
      lhXZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; //Gaussian kriging for intractable term
      
      logprob = -0.05 * thetaprop * thetaprop + 0.05 * thetaprev * thetaprev + lhXZp - lhXZ;              // log probability ratio to determine acceptance of MCMC 
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta(k+1) = thetaprop;
      lhXZ = lhXZp;		
    }else{
      theta(k+1) = thetaprev;
    }
    
  }
  
  return theta;
}




// [[Rcpp::export]]
// LikEm for Potts
vec pottsNormEm(int Niter, vec theta, double COV, double logconst, vec betahat, vec phihat, vec Designmat, vec y, double stat){
  int thnrow = Designmat.n_elem;                                        // number of design points
  double thetaprev, thetaprop;                                                   // befor propose in MCMC, previous parameters
  double logconstp,logprob,u;                                               // used in MCMC step
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], sigmasqhat = phihat[1]; 
  
  mat h1(thnrow,thnrow);   // used to construct Sigma in Gaussian process
  vec h1dcross(thnrow);     // used to construct cross Sigma in Gaussian process
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i)-Designmat(j));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat);
  mat InvSigma = inv(Sigma);	    // Inverse of Sigma in Gaussian process
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	// make design matrix for linear model 
  
  
  // Start of MCMC Chain 
  for(int k = 0; k< Niter-1; k++){
    
    double Znormal = randn(); // multivariate proposal by using Cholesky factorization
    thetaprev = theta(k);
    
    // proposed parameter and corresponding rho coefficients
    thetaprop = thetaprev + Znormal*COV;
    
    // constranits on prior space
    if( thetaprop < 0 ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  // Caculating cross covaraince matrix
        h1dcross[i] =  fabs(thetaprop-Designmat(i));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat);
      logconstp = (betahat(0) + betahat(1)*thetaprop + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; //Gaussian kriging for intractable term
      
      logprob = -0.05 * thetaprop * thetaprop + 0.05 * thetaprev * thetaprev + stat * (thetaprop - thetaprev) + logconst - logconstp;    // log probability ratio to determine acceptance of MCMC 
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta(k+1) = thetaprop;
      logconst = logconstp;			
    }else{
      theta(k+1) = thetaprev;
    }
    
  }
  
  return theta;
}




// =============================================================================
// Approximation to d
// =============================================================================

// [[Rcpp::export]]
// compute IS weights
vec compW(double thi, double thu, vec Sy){
  vec dummy = exp( Sy * (thi - thu) );
  return dummy / sum(dummy);
}



// [[Rcpp::export]]
mat cpp_appx(vec th, double Sx, double invcovth, vec thu, List aux){
  int nth = th.size();
  vec dummy, Sy, wIS;
  int indu;
  double Sybar, uhat, dhat;
  mat res = zeros(nth, 2);
  
  for(int thi = 0; thi < nth; thi ++){
    dummy = sqrt( (th[thi] - thu) % (th[thi] - thu) * invcovth );
    indu = dummy.index_min();
    Sy = as<vec>(aux[indu]);
    wIS = compW(th[thi], thu[indu], Sy);
    
    Sybar = sum( Sy % wIS );
    uhat = Sx - Sybar;
    dhat = -sum(wIS % Sy % Sy) + Sybar * Sybar + uhat * uhat;
    
    res(thi,0) = uhat;
    res(thi,1) = dhat;
  }
  
  return res;
}

