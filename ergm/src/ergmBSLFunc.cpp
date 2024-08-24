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




// [[Rcpp::export]]
double Choose(int x, int y){
  double result1 = 1, result2 = 1, result;
  int iter = y;
  
  if( x< y ){ result = 0;	}else{	
    for(int i = 0; i<iter; i++){
      result1 = result1*x;
      result2 = result2*y;
      y = y-1;
      x = x-1;   	
    }	
    
    
    result = result1/result2;
  }
  return(result);
}



// [[Rcpp::export]]
vec countShared(vec rowsum, mat X){
  
  int numedges = sum(rowsum)/2, nrow = rowsum.n_elem ,ind;
  vec hist = zeros(numedges);
  int num = 0;
  for(int k = 0; k< nrow; k++){
    for(int j = 0; j< (k+1); j++){ 
      if( X(k,j) == 1){  for(int i = 0; i<nrow; i++){ hist(num) =  hist(num) + X(i,k)*X(k,j)*X(j,i); }
      num = num + 1; }
      
    }
  }
  vec shared = zeros(nrow);
  for(int i = 0; i < numedges; i++ ){
    ind = hist(i);	
    if(ind>0){ shared(ind-1) = shared(ind-1) + 1; }
  }
  return(shared);
}



// [[Rcpp::export]]
vec Summary(mat X){
  int nrow = X.n_rows;
  double decay = 0.25;
  vec rowsum = sum(X,1), result = zeros(2);
  
  
  vec shared = countShared(rowsum,X);
  
  
  for(int i = 0; i< nrow; i++){
    result(0) = result(0) + Choose(rowsum(i),1);
    result(1) = result(1) + (   1- pow( (1-exp(-decay)),i+1)  )*shared(i);   
  }
  
  result(0)=result(0)/2;
  result(1)=exp(decay)*result(1);   
  
  return result;
}



// [[Rcpp::export]]
vec Gibbs(mat X, vec coef, int cycle){
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(2), Sumstat = Summary(X) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          
          changestat(1) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1);} 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          
          changestat(1) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej-1); 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        
        ivec(i) = 0; jvec(j) = 0;
      }
    }   
  }
  
  return(Sumstat); 
}


// [[Rcpp::export]]
mat Gibbs2(mat X, vec coef, int cycle){
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(2), Sumstat = Summary(X) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          
          changestat(1) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1);} 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          
          changestat(1) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej-1); 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        
        ivec(i) = 0; jvec(j) = 0;
      }
    }   
  }
  
  return(X); 
}


// [[Rcpp::export]]
mat Gibbs3(mat X, vec coef, int cycle){
  int nrow = X.n_rows,indi,indj,prei,prej,res;
  double decay = 0.25;
  vec star = sum(X,1), changestat = zeros(2), Sumstat = Summary(X) ;
  rowvec ivec = trans( zeros(nrow) ), jvec = trans( zeros(nrow) );
  vec geoweight = zeros(nrow);
  for(int i =0; i<nrow; i++){ geoweight(i) = 1- pow( (1-exp(-decay)),i+1); }
  
  mat result = zeros(cycle, 2);
  
  for(int l = 0; l< cycle; l++){
    for(int i = 1; i< nrow; i++){
      for(int j = 0; j< i; j++){
        ivec(i) = 1; jvec(j) = 1;
        res = 0; 
        
        
        if(X(i,j)==0){  
          indi = star(i) + 1, indj = star(j) + 1; 
          
          changestat(0) = ( Choose(indi,1) + Choose(indj,1) - Choose(star(i),1) - Choose(star(j),1) )/2;
          
          
          changestat(1) = 0;     
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei+1-1) ; 
              if(prei >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1);} 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej+1-1) ; 
              if(prej >0){ changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1);} 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() < p  ){
            X(i,j) = X(j,i) = 1; 
            star(i) = indi; star(j) = indj;		 
            Sumstat = Sumstat + changestat;
          }
          
          
        }else{
          indi = star(i) - 1, indj = star(j) - 1;	
          
          changestat(0) = ( Choose(star(i),1) + Choose(star(j),1) - Choose(indi,1) - Choose(indj,1) )/2;
          
          
          changestat(1) = 0;
          for(int k = 0; k< nrow; k++){
            if(   ( X(i,k) == 1 ) & ( X(j,k) == 1  )   ){
              prei = sum( X.row(i)%X.row(k) );
              prej = sum( X.row(j)%X.row(k) );
              if(prei-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prei-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prei-1); 
              if(prej-1 >0){changestat(1) = changestat(1) - exp(decay)*geoweight(prej-1-1); } 
              changestat(1) = changestat(1) + exp(decay)*geoweight(prej-1); 
              res = res + 1; 
            }
          }
          if(res > 0){ changestat(1) = changestat(1) + exp(decay)*geoweight(res-1); }
          
          
          vec r =  exp(trans(coef)*changestat) ;         
          double p = r(0)/(1+r(0));   		   
          if( randu() > p  ){
            X(i,j) = X(j,i) = 0; 
            star(i) = indi; star(j) = indj;  
            Sumstat = Sumstat - changestat;
          }             
        }
        
        ivec(i) = 0; jvec(j) = 0;
      }
    }
    result(l,0) = Sumstat[0];
    result(l,1) = Sumstat[1];
  }
  
  return(result); 
}


// [[Rcpp::export]]
cube pGibbs3(mat X, vec coef, int cycle, int m, int num){
  
  cube H(cycle,m,2);                       
  omp_set_num_threads(num);
  
#pragma omp parallel num_threads(num)
{
#pragma omp for
  for(int M = 0; M < m; M ++){
    mat dummy = Gibbs3(X, coef, cycle);
    
    for(int l = 0; l < cycle; l ++){
      H(l,M,0) = dummy(l,0);
      H(l,M,1) = dummy(l,1);
    }
  }
}
return(H);
}


// [[Rcpp::export]]
cube pAuxSamp(mat X, int cycle, mat Designmat, int m, int num){
  
  int thnrow = Designmat.n_rows;            
  cube H(thnrow,2,m);                       
  omp_set_num_threads(num);
  
  for(int M = 0; M < m; M++){            
    
    int i;
#pragma omp parallel shared(Designmat) private(i)
{	
#pragma omp for schedule(static)  
  for(i = 0; i < thnrow; i++){
    vec thetaprop = trans( Designmat.row( i ) );        
    vec sumstat = Gibbs(X, thetaprop, cycle);
    
    for(int j = 0; j < 2; j++){ H(i,j,M) = sumstat[j];  }
    
  }
}
  }
  return(H);        	
}


// [[Rcpp::export]]
mat pResponseErgm(mat X, int cycle, vec hatparameter, int m, int num){
  mat H(m,2);                       
  omp_set_num_threads(num);
  
  int M;
#pragma omp parallel shared(H) private(M)
{	
#pragma omp for schedule(static)  
  for(M = 0; M < m; M++){                           
    vec sumstat = Gibbs(X, hatparameter, cycle);  
    H(M,0) = sumstat[0];
    H(M,1) = sumstat[1];
  }
}
return(H);        	
}


// [[Rcpp::export]]
mat ergmDMH(mat X, mat COV, mat theta, int outer, int cycle){
  
  double logprob,u;                        
  int nCOVcols = COV.n_cols;                
  vec thetaprev(nCOVcols);                  
  vec stat = Summary(X), statprop(nCOVcols); 
  
  
  for(int l = 0; l< outer; l++){
    
    if( (l > 1000) && (l <= 10000) ){ 
      COV = cov(theta);
    }	
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  
    
    
    statprop = Gibbs(X, thetaprop, cycle);
    
    
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + trans(thetaprev - thetaprop)*(statprop - stat) );
    logprob = dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
  }
  
  return(theta);	
}


// [[Rcpp::export]]
mat ergmFDMH(mat X, mat COV, mat theta, int outer, int cycle){
  
  double logprob,u;                        
  int nCOVcols = COV.n_cols;                
  vec thetaprev(nCOVcols);                  
  vec stat = Summary(X), statprop(nCOVcols); 
  
  
  for(int l = 0; l< outer; l++){
    
    if( (l > 1000) && (l <= 10000) ){ 
      COV = cov(theta);
    }	
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  
    
    
    statprop = Gibbs(X, thetaprop, cycle);
    
    
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + trans(thetaprev - thetaprop)*(statprop - stat) );
    logprob = 0.5 * dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
  }
  
  return(theta);	
}


// [[Rcpp::export]]
mat ergmAEX(int Niter, int Numaux, int cycle, double t0, int neighbor, mat th, mat thdist, mat theta, mat COV, mat X){ 
  int auxNiter = Numaux*20 + Numaux;
  double d = th.n_rows;                                        // number of particles
  vec p = (1/d)*ones(d), lw = zeros(d), Vis = zeros(d);        // target prob, logartihm of abundance factor
  mat  Sdummy = zeros(auxNiter,5);                             // data base suff, parameter, abundance factor will be stored
  
  // initialize
  int indprop, ind = 0;
  Sdummy(0,4) = lw[ind];
  Sdummy(0, span(2,3)) = th.row( ind );
  mat auxvar = Gibbs2(X, trans(th.row( ind )), cycle);
  Sdummy(0, span(0,1)) = trans(  Summary( auxvar )  );
  
  // preliminary run of the auxiliary chain 
  for(int k = 0; k< auxNiter-1; k++){
    // decide update theta or aux var 
    if( randu() < 0.75 ){  // update theta
      uvec dummy = sort_index(  thdist.row( ind )   );  // using distance matrix, ordering from current theta 
      uvec samplist = dummy( span( 1, neighbor )  );    // neighborhood calcul
      int ii  = floor( (neighbor)*randu() ) ;         // sample uniformly from neighborhood  
      indprop = samplist[ii];
      
      rowvec thprop =  th.row( indprop );
      Sdummy(k+1, span(0,1)) = Sdummy(k, span(0,1));
      vec caldummy = Sdummy(k, span(0,1))*trans( thprop - Sdummy(k, span(2,3)) );
      
      double logprob = lw[ind] - lw[indprop]  + ( caldummy )[0]; 
      double u = log( randu() );
      if( u< logprob ){
        Sdummy(k+1, span(2,3)) = thprop;
        ind = indprop;
      }else{
        Sdummy(k+1, span(2,3)) = Sdummy(k, span(2,3));
        ind = ind;  
      }  
      
    }else{                // update aux var, always accept by DBE
      auxvar = Gibbs2(auxvar, trans(th.row( ind )), cycle);  
      ind = ind, Sdummy(k+1, span(2,3)) = Sdummy(k, span(2,3)), Sdummy(k+1, span(0,1)) = trans( Summary(auxvar) );
    }
    
    // update abundance factor
    vec e = zeros(d);
    e[ind] = 1;
    double val;
    double kdummy = k;
    if( t0 > kdummy ){ val = t0; }else{ val = kdummy; }
    lw = lw + (t0/val)*(e-p);
    Sdummy(k+1,4) = lw[ind];
    
    // record visitation frequencies
    Vis = Vis + e;
  }  
  // Rprintf("Finished preliminary run of auxiliary chain...");
  
  
  // burn in Numaux number and equally 20 sample spaced
  mat S = zeros(Numaux+Niter,5); 
  for(int k = 0; k< Numaux; k++){ S.row(k) = Sdummy.row(Numaux-1 + (k+1)*20); }
  
  // initialize final chain
  double logprob;                             // used in Outer MCMC
  int nCOVcols = COV.n_cols;                 // number of parameters 
  vec thetaprev(nCOVcols);                   // before propose in Outer MCMC, previous parameters
  vec XStat = Summary(X);                    // sufficient statistics
  mat parameter(Niter,nCOVcols);             // posterior samples will be stored          
  parameter.row( 0 ) = theta.row(0);
  
  // run auxiliary chain and target chain simultaneously
  for(int k = 0; k< Niter-1; k++){
    
    // if( k >= 200000-1 ){ Rprintf("%dth iteration...\n", k+1); }
    // Rprintf("%dth iteration...\n", k+1);
    
    // Auxiliary chain 
    // decide update theta or aux var 
    if( randu() < 0.75 ){  // update theta
      uvec dummy = sort_index(  thdist.row( ind )   );  // using distance matrix, ordering from current theta 
      uvec samplist = dummy( span( 1, neighbor )  );    // neighborhood calcul
      int ii  = floor( (neighbor)*randu() ) ;         // sample uniformly from neighborhood  
      indprop = samplist[ii];
      
      rowvec thprop =  th.row( indprop );
      S(Numaux-1 + k + 1, span(0,1)) = S(Numaux-1 + k, span(0,1));
      vec caldummy = S(Numaux-1 + k, span(0,1))*trans( thprop - S(Numaux-1 + k, span(2,3)) );
      
      double logprob1 = lw[ind] - lw[indprop]  + ( caldummy )[0]; 
      double u = log( randu() );
      if( u< logprob1 ){
        S(Numaux-1 + k + 1, span(2,3)) = thprop;
        ind = indprop;
      }else{
        S(Numaux-1 + k + 1, span(2,3)) = S(Numaux-1 + k, span(2,3));
        ind = ind;  
      }                   
    }else{                // update aux var, always accept by DBE
      auxvar = Gibbs2(auxvar, trans(th.row( ind )), cycle);  
      ind = ind, S(Numaux-1 + k + 1, span(2,3)) = S(Numaux-1 + k, span(2,3)), S(Numaux-1 + k + 1, span(0,1)) = trans( Summary(auxvar) );
    }
    // if( k >= 200000 ){ Rprintf("Finished auxilary chain...\n"); }
    // Rprintf("Finished auxilary chain...\n"); 
    
    // update abundance factor
    vec e = zeros(d);
    e[ind] = 1;
    double val;
    double kdummy = auxNiter-1 + k;
    if( t0 > kdummy ){ val = t0; }else{ val = kdummy; }
    lw = lw + (t0/val)*(e-p);
    S(Numaux-1 + k + 1, 4) = lw[ind];
    // if( k >= 200000 ){ Rprintf("Finished updating abundance factor...\n"); }
    // Rprintf("Finished updating abundance factor...\n"); 
    
    
    // Target chain
    if( (k > 1000) && (k <= 10000) ){ // adaptively update COV until 10000 iterations 
      COV = cov(  parameter(span(0,k),span(0,1))  );
    } 
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = parameter(k,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           // multivariate proposal by using Cholesky factorization
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  // proposed parameter
    
    vec Prob = zeros(Numaux + k + 1);  // because the number of row of S is Numaux + k + 1
    for(int m = 0; m< Numaux + k + 1; m++){ 
      vec caldummy2 = S(m,4) +  S(m, span(0,1))*(   thetaprop - trans( S(m, span(2,3)) )   );   
      Prob[m] = caldummy2[0]  ;
    }
    
    double mx = (Prob).max();
    Prob = exp(Prob-mx);
    Prob = Prob/sum(Prob);     
    double uu = randu(), l = -1, om = 0;
    while( om<uu ){
      l = l+1;
      om = om + Prob[l];
    }
    vec Statprop = trans( S(l,span(0,1)) );  
    vec caldummy3 = trans(thetaprop - thetaprev)*XStat + trans(thetaprev - thetaprop)*Statprop;
    
    logprob = caldummy3[0];   
    if( log( randu() )< logprob ){
      parameter.row( k+1 ) = trans(thetaprop); 
    }else{
      parameter.row( k+1 ) = trans(thetaprev); 
    }
    // if( k >= 200000 ){ Rprintf("Finished target chain...\n"); }
    // Rprintf("Finished target chain...\n");
  }
  
  // return List::create(Named("Vis") = Vis, 
  //                     Named("S") = S, 
  //                     Named("abundance") = lw, 
  //                     Named("par") = parameter); 
  
  return parameter;
}


// [[Rcpp::export]]
//  Youne's stochastic approximation #
mat SApprox(mat X, mat Domain, mat th, int cycle, int Niter){
  
  double gam0 = 1, eps = 0.1;
  vec Stat0 = Summary(X), Stat;
  mat Data = X;
  int d = th.n_rows, p = th.n_cols;        
  
  for(int l = 0; l< d; l++){
    vec theta = trans(th.row( l ));
    X = Data;
    Stat = Stat0;
    
    for(int i = 0; i< Niter; i++){
      double gam = gam0;
      Stat = Gibbs(X, theta, cycle);
      theta = theta + gam*(Stat0-Stat);
      
      for(int j = 0; j<p; j++){
        if(theta[j]<Domain(0,j)){ theta[j] = Domain(0,j) + eps; }
        if(theta[j]>Domain(1,j)){ theta[j] = Domain(1,j) - eps; }	
      }			
    }
    th.row( l ) =  trans(theta);		
  }	  
  
  return(th);
}


// [[Rcpp::export]]
// Atchade's adpative algorithm for ergm
mat ergmAtchade(int outer, int cycle, mat COV, mat theta, mat th, mat X){	    
  int d = th.n_rows, p = th.n_cols, Vmin;  
  vec Vis= zeros(d), gamVec(4), cVec;    // gamVec length and sequence both can be changed 
  for(int i = 0; i< 4; i++){ gamVec[i] = pow( 0.1, i ); }
  
  mat Data = X;
  vec Stat0 = Summary(X), Stat;
  // summary statistics will be stored
  mat E1sum = zeros(outer,d), E2sum = zeros(outer,d); 
  
  // approximate cVec ( logZ(theta_{i}) ) until gam is small
  for(int igam = 0; igam< 4; igam++){  
    
    if( igam == 3 ){ Vmin=1000; }else{ Vmin=0; } //Vmin=1000 can be changed
    double gam = gamVec(igam);
    X = Data, Stat = Stat0; 
    vec VisTemp = zeros(d);     
    
    // pos denotes the variable I in Algorithm 2.1
    int pos = d-1;                          
    
    // cVec is the variable c in Algorithm 2.1
    if( igam == 0 ){ cVec = zeros(d); }
    
    // Stopping criteria
    while( VisTemp.min() <= Vmin || abs(  VisTemp-mean(VisTemp)*ones(d)  ).max() >= 0.2*mean(VisTemp) ){        
      
      // Update X_{n}
      Stat = Gibbs(Data, trans(th.row( pos )), cycle);
      
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
  vec thetaprev(p), thetaprop(p);
  X = Data, Stat = Stat0;
  int pos = d-1;    
  double gam = gamVec(3),prob; 
  
  // c = cVec from above iteration
  // Bandwidth for smoothing the estimate
  int hpar = 20;
  for(int k = 0; k< outer-1; k++){
    if( (k > 1000) && (k <= 10000) ){ // adaptively update COV until 10000 iterations 
      COV = cov(theta);
    }	
    
    for(int i = 0; i< p; i++){
      thetaprev[i] = theta(k,i);
    }
    // proposed parameter   
    vec Znormal = randn(p);                                                  
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV + 0.001 * diagmat(ones(p))) ); 
    
    // Update X_{n}
    Stat = Gibbs(Data, trans(th.row( pos )), cycle);
    
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
    Vis[pos] = Vis[pos] + 1;
    E1sum(Vis[pos]-1,pos) = Stat[0], E2sum(Vis[pos]-1,pos) = Stat[1];
    
    // Update theta if there are enough observations
    // Evaluate \pi(theta)
    // Distance to each theta^i
    mat dummy = thetaprev[0]*ones(d);
    dummy.insert_cols( 1, thetaprev[1]*ones(d) );
    uvec ind = sort_index(  sum( pow( (dummy - th), 2 ), 1)  );
    vec cw = (1/hpar)*ones(hpar);
    
    // There might be some trouble here if none of the closest 'hpar'
    // particles around thetaprop has received data
    vec expEtheta = zeros(hpar);
    for(int i = 0; i< hpar; i++){	
      if( Vis[ind[i]] == 0 ){ expEtheta[i] = 0; }else{
        vec tmp = ( thetaprev[0] - th(ind[i],0) )*E1sum.col( ind[i] ) + ( thetaprev[1] - th(ind[i],1) )*E2sum.col( ind[i] );	
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
    double Eth = Stat0[0]*thetaprev[0] + Stat0[1]*thetaprev[1] - ftheta - log(sum( exp(dummy2 - ftheta) ));
    
    
    // Evaluate posterior at the new theta        
    dummy = thetaprop[0]*ones(d);
    dummy.insert_cols( 1, thetaprop[1]*ones(d) );
    ind = sort_index(  sum( pow( (dummy - th), 2 ), 1)  );        
    cw = (1/hpar)*ones(hpar);
    // There might be some trouble here if none of the closest 'hpar'
    // particles around thetaprop has received data
    vec expEthetaprop = zeros(hpar);        
    for(int i = 0; i< hpar; i++){	
      if( Vis[ind[i]] == 0 ){ expEthetaprop[i] = 0; }else{
        vec tmp = ( thetaprop[0] - th(ind[i],0) )*E1sum.col( ind[i] ) + ( thetaprop[1] - th(ind[i],1) )*E2sum.col( ind[i] );  
        tmp = tmp(  span( 0, Vis[ind[i]]-1 )  );
        double tmpMax = tmp.max();	
        expEthetaprop[i] = tmpMax + log( sum( exp( tmp - tmpMax ) ) ) - log( Vis[ind[i]] ); // eq(8)'s [   ] part       	
      }
    }
    for(int i = 0; i< hpar; i++){ dummy2[i] = cVec[ ind( span(0,hpar-1) )[i] ] + expEthetaprop[i] + cw[i];	}
    double fthetaprop = dummy2.max();   
    double Ethprop = Stat0[0]*thetaprop[0] + Stat0[1]*thetaprop[1] - fthetaprop - log(sum( exp(dummy2 - fthetaprop) ));						
    
    // Acceptance prob.
    prob = exp(Ethprop-Eth);
    
    // Accept - Reject
    u = randu();
    if(u <= prob){ 
      theta.insert_rows(k+1,trans(thetaprop));
    }else{
      theta.insert_rows(k+1,trans(thetaprev));	
    }
  }
  
  return(theta);
}



// [[Rcpp::export]]
mat GPmcmcErgmLik(int Niter, mat theta, mat COV, double lhXZ, vec betahat, vec phihat, mat Designmat, vec y, vec stat){
  int thnrow = Designmat.n_rows;                                                
  int nCOVcols = COV.n_cols;                                                    
  vec thetaprev(nCOVcols);                                                     
  
  double lhXZp = 0, logprob = 0, u = 0;                                               
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], phi2hat = phihat[1], sigmasqhat = phihat[2]; 
  
  
  int percentile = 0.0025*thnrow; 
  
  
  mat Domain(2,nCOVcols);                                                         
  for(int i = 0; i < nCOVcols; i++){
    vec dummy = sort( Designmat.col( i ) );
    Domain(0,i) = dummy(percentile);
    Domain(1,i) = dummy(thnrow-1-percentile);
  }
  
  mat h1(thnrow,thnrow), h2(thnrow,thnrow);   
  vec h1dcross(thnrow), h2dcross(thnrow);     
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i,0)-Designmat(j,0));
      h2(i,j) = h2(j,i) = fabs(Designmat(i,1)-Designmat(j,1));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat)%(1+sqrt(3)*h2/phi2hat)%exp(-sqrt(3)*h2/phi2hat);
  mat InvSigma = inv(Sigma);	   
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	 
  
  
  
  for(int k = 0; k< Niter; k++){
    if( (k > 1000) && (k <= 10000) ){ 
      COV = cov(theta);
    }	
    vec Znormal = randn(nCOVcols);
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(k,i);
    }
    
    
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );
    
    
    if( thetaprop[0] > Domain(1,0) || thetaprop[0] < Domain(0,0) || thetaprop[1] > Domain(1,1) || thetaprop[1] < Domain(0,1) ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  
        h1dcross[i] =  fabs(thetaprop[0]-Designmat(i,0));
        h2dcross[i] =  fabs(thetaprop[1]-Designmat(i,1));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat)%(1+sqrt(3)*h2dcross/phi2hat)%exp(-sqrt(3)*h2dcross/phi2hat);
      vec xpoint = ones(1);
      xpoint.insert_rows(1,thetaprop);
      lhXZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; 
      
      logprob = lhXZp - lhXZ;             
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(k+1,trans(thetaprop));
      lhXZ = lhXZp;		
    }else{
      theta.insert_rows(k+1,trans(thetaprev));
    }
    
  }
  
  return theta;
}


// [[Rcpp::export]]
mat GPmcmcErgmLikPrior(int Niter, mat theta, mat COV, double lhXZ, vec betahat, vec phihat, mat Designmat, vec y, vec stat){
  int thnrow = Designmat.n_rows;                                                
  int nCOVcols = COV.n_cols;                                                    
  vec thetaprev(nCOVcols);                                                     
  
  double lhXZp = 0, logprob = 0, u = 0;                                               
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], phi2hat = phihat[1], sigmasqhat = phihat[2]; 
  
  
  // int percentile = 0.0025*thnrow; 
  // 
  // 
  // mat Domain(2,nCOVcols);                                                         
  // for(int i = 0; i < nCOVcols; i++){
  //   vec dummy = sort( Designmat.col( i ) );
  //   Domain(0,i) = dummy(percentile);
  //   Domain(1,i) = dummy(thnrow-1-percentile);
  // }
  
  mat h1(thnrow,thnrow), h2(thnrow,thnrow);   
  vec h1dcross(thnrow), h2dcross(thnrow);     
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i,0)-Designmat(j,0));
      h2(i,j) = h2(j,i) = fabs(Designmat(i,1)-Designmat(j,1));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat)%(1+sqrt(3)*h2/phi2hat)%exp(-sqrt(3)*h2/phi2hat);
  mat InvSigma = inv(Sigma);	   
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	 
  
  
  
  for(int k = 0; k< Niter; k++){
    if( (k > 1000) && (k <= 10000) ){ 
      COV = cov(theta);
    }	
    vec Znormal = randn(nCOVcols);
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(k,i);
    }
    
    
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );
    
    
    if( thetaprop[0] > 2.27 || thetaprop[0] < -5 || thetaprop[1] > 2.32 || thetaprop[1] < -1.57 ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  
        h1dcross[i] =  fabs(thetaprop[0]-Designmat(i,0));
        h2dcross[i] =  fabs(thetaprop[1]-Designmat(i,1));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat)%(1+sqrt(3)*h2dcross/phi2hat)%exp(-sqrt(3)*h2dcross/phi2hat);
      vec xpoint = ones(1);
      xpoint.insert_rows(1,thetaprop);
      lhXZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; 
      
      logprob = lhXZp - lhXZ;             
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(k+1,trans(thetaprop));
      lhXZ = lhXZp;		
    }else{
      theta.insert_rows(k+1,trans(thetaprev));
    }
    
  }
  
  return theta;
}



// [[Rcpp::export]]
mat GPmcmcErgmNorm(int Niter, mat theta, mat COV, double lZ, vec betahat, vec phihat, mat Designmat, vec y, vec stat){
  int thnrow = Designmat.n_rows;                                                
  int nCOVcols = COV.n_cols;                                                    
  vec thetaprev(nCOVcols);                                                     
  
  double lZp = 0, logprob = 0, u = 0;                                               
  double negativeInf = -std::numeric_limits<float>::infinity();;	               
  double phi1hat = phihat[0], phi2hat = phihat[1], sigmasqhat = phihat[2]; 
  
  int percentile = 0.0025*thnrow; 
  
  mat Domain(2,nCOVcols);                                                         
  for(int i = 0; i < nCOVcols; i++){
    vec dummy = sort( Designmat.col( i ) );
    Domain(0,i) = dummy(percentile);
    Domain(1,i) = dummy(thnrow-1-percentile);
  }
  
  mat h1(thnrow,thnrow), h2(thnrow,thnrow);   
  vec h1dcross(thnrow), h2dcross(thnrow);     
  
  for(int i = 0; i < thnrow; i++){
    for(int j = 0; j <= i; j++){
      h1(i,j) = h1(j,i) = fabs(Designmat(i,0)-Designmat(j,0));
      h2(i,j) = h2(j,i) = fabs(Designmat(i,1)-Designmat(j,1));
    }
  }
  mat Sigma = sigmasqhat*(1+sqrt(3)*h1/phi1hat)%exp(-sqrt(3)*h1/phi1hat)%(1+sqrt(3)*h2/phi2hat)%exp(-sqrt(3)*h2/phi2hat);
  mat InvSigma = inv(Sigma);	   
  mat Xth = ones(thnrow,1);
  Xth.insert_cols(1,Designmat);	 
  
  
  
  for(int k = 0; k< Niter; k++){
    
    if( (k > 1000) && (k <= 10000) ){ 
      COV = cov(theta);
    }	
    vec Znormal = randn(nCOVcols);
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(k,i);
    }
    
    
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );
    
    
    if( thetaprop[0] > Domain(1,0) || thetaprop[0] < Domain(0,0) || thetaprop[1] > Domain(1,1) || thetaprop[1] < Domain(0,1) ){
      logprob = negativeInf;	
      
    }else{			
      for(int i = 0; i< thnrow; i++){  
        h1dcross[i] =  fabs(thetaprop[0]-Designmat(i,0));
        h2dcross[i] =  fabs(thetaprop[1]-Designmat(i,1));	
      }
      mat Sigmacross = sigmasqhat*(1+sqrt(3)*h1dcross/phi1hat)%exp(-sqrt(3)*h1dcross/phi1hat)%(1+sqrt(3)*h2dcross/phi2hat)%exp(-sqrt(3)*h2dcross/phi2hat);
      vec xpoint = ones(1);
      xpoint.insert_rows(1,thetaprop);
      lZp = (trans(xpoint)*betahat + trans(Sigmacross)* InvSigma*(y-Xth*betahat))[0]; 
      
      logprob = (trans(thetaprop - thetaprev)*stat + (lZ - lZp))[0];
    } 
    
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(k+1,trans(thetaprop));
      lZ = lZp;		
    }else{
      theta.insert_rows(k+1,trans(thetaprev));
    }
    
  }
  
  return theta;
}



// [[Rcpp::export]]
// compute IS weights
vec compW(vec thi, vec thu, mat Sy){
  vec dummy = exp( Sy * (thi - thu) );
  return dummy / sum(dummy);
}



// [[Rcpp::export]]
mat rcppf_Hhat(vec wIS, mat Sy, vec theta){
  int p = theta.size(); 
  mat Hhat(p, p);
  
  for(int j = 0; j < p; j ++){
    for(int i = 0; i < p; i ++){
      if( i <= j ){
        Hhat(i,j) = - sum( wIS % Sy.col(i) % Sy.col(j) ) + sum( wIS % Sy.col(i) ) * sum( wIS % Sy.col(j) );
        if(i < j){
          Hhat(j,i) = Hhat(i,j);  
        }
      } 
    }  
  }
  
  return Hhat;
}





// [[Rcpp::export]]
mat cpp_appx(mat th, vec Sx, mat invcovth, mat thu, List aux){
  int nth = th.n_rows, p = th.n_cols, m = thu.n_rows, dummyind;
  vec dummy(m), wIS, uhat;
  mat Sy;
  int indu;
  mat res = zeros(nth, p + p*(p+1)/2);
  
  for(int thi = 0; thi < nth; thi ++){
    
    for(int j = 0; j < m; j ++){
      dummy[j] = sqrt( (th.row(thi) - thu.row(j)) * invcovth * trans(th.row(thi) - thu.row(j)) )[0];
    }
    indu = dummy.index_min();
    Sy = as<mat>(aux[indu]);
    
    wIS = compW(trans(th.row(thi)), trans(thu.row(indu)), Sy);
    uhat = Sx - trans(Sy) * wIS;
    res(thi, span(0, p-1)) = trans(uhat);
    
    dummyind = 0;
    for(int k = 0; k < p; k ++){
      for(int l = 0; l < p; l ++){
        if(k <= l){
          res(thi, p+dummyind) = - sum( wIS % Sy.col(k) % Sy.col(l) ) + 
            sum( wIS % Sy.col(k) ) * sum( wIS % Sy.col(l) ) +
            uhat[k] * uhat[l];
          dummyind = dummyind + 1;
        }
      }
    }
  }
  
  return res;
}




// [[Rcpp::export]]
vec AIKS(mat th, mat score, vec weight, double c, double beta, int i, int k){
  int p = score.n_cols;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip;
  mat knot = zeros(k,p);
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    }
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    }
  }
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    double wsq = pow(weight[i],2)*knot(i,j);
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        wsq += weight[i]*knot(m,j)*weight[m]*2;
      }
    }
    H[j] = wsq;
  }
  
  return(H);
}



// [[Rcpp::export]]
// Compute w^2 for each row
mat pAIKS(mat th, mat score, vec weight, double c, double beta, int k, int num){
  int p = score.n_cols;
  mat H(k,p);
  omp_set_num_threads(num);
  
  int i;
#pragma omp parallel shared(H) private(i)
{
#pragma omp for schedule(static)
  for(i = 0; i < k; i++){
    vec res = AIKS(th, score, weight, c, beta, i, k);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
}
return(H);
}



// [[Rcpp::export]]
// Compute w^2 for each row
mat npAIKS(mat th, mat score, vec weight, double c, double beta, int k){
  int p = score.n_cols;
  mat H(k,p);
  
  for(int i = 0; i < k; i++){
    vec res = AIKS(th, score, weight, c, beta, i, k);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
  return(H);
}




// [[Rcpp::export]]
// For bootstrapping
vec AIKSboot(mat th, mat score, vec rep, vec auxdep, double c, double beta, int i, int k){
  int p = score.n_cols, nn = rep.size(), nind, nind2;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip, wsq;
  mat knot = zeros(k,p);
  uvec ind, ind2;
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    }
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    }
  }
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    
    ind = find(rep == i);
    nind = ind.size();
    
    wsq = 0;
    for(int m1 = 0; m1< nind; m1 ++){
      for(int m2 = 0; m2 < nind; m2 ++){
        wsq += auxdep[ ind[m1] ] * auxdep[ ind[m2] ] * knot(i,j) / pow(nn, 2);
      }
    }
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        
        ind2 = find(rep == m);
        nind2 = ind2.size();
        
        for(int m1 = 0; m1< nind; m1 ++){
          for(int m2 = 0; m2 < nind2; m2 ++){
            wsq +=  auxdep[ ind[m1] ] * auxdep[ ind2[m2] ] * knot(m,j) * 2 / pow(nn, 2);
          }
        }  
      }
    }
    H[j] = wsq;
  }
  
  return(H);
}





// [[Rcpp::export]]
// Compute w^2 for each row
mat npAIKSboot(mat th, mat score, vec rep, vec auxdep, double c, double beta, int k){
  int p = score.n_cols;
  mat H(k,p);
  
  for(int i = 0; i < k; i++){
    vec res = AIKSboot(th, score, rep, auxdep, c, beta, i, k);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
  return(H);
}



// [[Rcpp::export]]
vec AIKS_woWeight(mat th, mat score, double c, double beta, int i, int k, int num){
  int p = score.n_cols;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip;
  mat knot = zeros(k,p);
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    }
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    }
  }
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    double wsq = knot(i,j);
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        wsq += knot(m,j)*2;
      }
    }
    H[j] = wsq;
  }
  
  return(H);
}





// [[Rcpp::export]]
// Compute w^2 for each row
mat pAIKS_woWeight(mat th, mat score, double c, double beta, int k, int num){
  int p = score.n_cols;
  mat H(k,p);
  omp_set_num_threads(num);
  
  int i;
#pragma omp parallel shared(H) private(i)
{
#pragma omp for schedule(static)
  for(i = 0; i < k; i++){
    vec res = AIKS_woWeight(th, score, c, beta, i, k, num);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
}
return(H);
}




// [[Rcpp::export]]
mat ergmBSL(mat X, mat COV, mat theta, int outer, int auxiter, int cycle){
  
  double logprob,u;                        
  int nCOVcols = COV.n_cols;                
  vec thetaprev(nCOVcols);
  vec stat = Summary(X); 
  mat auxprop(auxiter, nCOVcols); //samples of summary statistics
  mat auxprev(auxiter, nCOVcols); //samples of summary statistics
  vec propstat(nCOVcols);
  vec prevstat(nCOVcols); // auxiliary sample
  
  
  for(int l = 0; l< outer; l++){
    
    if( (l > 1000) && (l <= 10000) ){
      COV = cov(theta);
    }
    
    for(int i = 0; i< nCOVcols; i++){
      thetaprev[i] = theta(l,i);
    }
    
    vec Znormal = randn(nCOVcols);                                           
    vec thetaprop = trans(  trans(thetaprev) + trans(Znormal)*chol(COV)  );  
    
    for(int m = 0; m < auxiter; m++){
      // sample with proposal
      propstat = Gibbs(X, thetaprop, cycle);
      auxprop.row(m) = trans(propstat);
      
      // sample with previous
      prevstat = Gibbs(X, thetaprev, cycle);
      auxprev.row(m) = trans(prevstat);
    }
    
    vec mprop = trans(mean(auxprop, 0));
    vec mprev = trans(mean(auxprev, 0));
    mat Vprop = cov(auxprop);
    mat Vprev = cov(auxprev);

    // Regularization
    double regularization = 1e-6;

    if(!Vprop.is_sympd()){ //if not positive definite
      Vprop = Vprop + regularization * eye(size(Vprop));
    }
    if(!Vprev.is_sympd()){
      Vprev = Vprev + regularization * eye(size(Vprev));;
    }

    mat InvVprop = inv(Vprop);
    mat InvVprev = inv(Vprev);
    
    vec llprop = -0.5*log(det(Vprop)) - 0.5*trans(stat-mprop)*InvVprop*(stat-mprop);
    vec llprev = -0.5*log(det(Vprev)) - 0.5*trans(stat-mprev)*InvVprev*(stat-mprev);
    
    vec dummy = ( -0.05*trans(thetaprop)*thetaprop + 0.05*trans(thetaprev)*thetaprev + llprop - llprev);
    logprob = dummy[0];
    u = log( randu() );
    if( u< logprob ){
      theta.insert_rows(l+1,trans(thetaprop));
    }else{
      theta.insert_rows(l+1,trans(thetaprev));
    }
    
  }
  
  return(theta);	
}



