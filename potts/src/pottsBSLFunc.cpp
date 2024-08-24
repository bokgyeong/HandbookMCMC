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



// [[Rcpp::depends("RcppArmadillo")]]
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



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generate summary statistics of auxiliary variables for given particles to construct IS estimate

vec pResponse(RawVector foo, int ncolor,  int cycle, double hatparameter, int m, int num){

vec H(m);                         // matrix where summary statistics will be stored. (m by 2)
omp_set_num_threads(num);

//////    START OF BIGGEST CHAIN (M)     ////// 
int M;
#pragma omp parallel shared(H) private(M)
{	
#pragma omp for schedule(static)  
for(M = 0; M < m; M++){                            // m is number of importance sampling estimate 
    // summary statistics
    H(M) = potts_stat(foo, ncolor, hatparameter, cycle);
}
}

return(H);        	
}



// [[Rcpp::depends("RcppArmadillo")]]
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



// [[Rcpp::depends("RcppArmadillo")]]
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



// [[Rcpp::export]]
List pottsBSL(RawVector foo, int ncolor, double stat, double COV, double beta, int outer, int auxiter, int cycle,
              bool updateCOV, double sigma2, int adaptInterval, double adaptFactorExponent, int adapIter){
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double logprob, u, betaprev, betaprop, statprop, statprev;
  vec postSamples = zeros(outer), accprob = zeros(outer);
  double rhat, gamma1, gamma2, c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double cholCOV = sqrt( sigma2 * COV );;
  vec auxprop(auxiter); //samples of summary statistics
  vec auxprev(auxiter); //samples of summary statistics
  
  for(int l = 0; l< outer; l++){
    
    //adaptive update of stepsize
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
    
    // propose new beta
    betaprev = beta;
    betaprop = betaprev + cholCOV * randn();
    
    if(betaprop <= 0){
      logprob = negativeInf;
    }else{
      // make auxiliary samples
      for(int m = 0; m < auxiter; m++){
        // sample with proposal
        statprop = potts_stat(foo, ncolor, betaprop, cycle);
        auxprop(m) = statprop;
        
        // sample with previous
        statprev = potts_stat(foo, ncolor, betaprev, cycle);
        auxprev(m) = statprev;
      }
      
      // construct S MVN(sample mean, sample COV)
      double mprop = mean(auxprop);
      double Vprop = var(auxprop);
      double mprev = mean(auxprev);
      double Vprev = var(auxprev);
      
      double llprop = -0.5*log(Vprop) - 0.5*(1/Vprop)*pow(stat-mprop,2);
      double llprev = -0.5*log(Vprev) - 0.5*(1/Vprev)*pow(stat-mprev,2);
      
      // accept or reject
      logprob = -0.05*betaprop*betaprop + 0.05*betaprev*betaprev + llprop - llprev;
    }
      u = log( randu() );
      if( u < logprob ){
        beta = betaprop;
        accprob[l] = 1;
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
