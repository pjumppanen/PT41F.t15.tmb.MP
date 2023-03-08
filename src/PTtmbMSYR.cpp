// Pella Tomlinson production model
// derived from 

// Space time
#include <TMB.hpp>

// square
template<class Type>
Type square(Type x){
  return pow(x,2.0); 
}

// sqrt
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5); 
}



// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( I_t ); // CPUE               
  DATA_VECTOR( c_t );                
  DATA_VECTOR( priorMode );
  DATA_VECTOR( priorLogCV );
  DATA_VECTOR( log_kBounds );
  
  // Parameters
  PARAMETER( log_MSY );
  PARAMETER( log_k );
  //PARAMETER( log_q );
  PARAMETER( shape );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaI );
  PARAMETER_VECTOR( log_B_t );  //Biomass with the productivity deviation
  
  // Derived quantities
  int n_years = I_t.size();
  
  Type k = exp(log_k); //unbounded k
  Type r = exp(log_MSY)*pow((shape+1),(1/shape))/k;

  //Type q = exp(log_q);
  vector<Type> Bpred_t( n_years ); //Biomass without the productivity deviation
  vector<Type> recDev( n_years ); //Biomass without the productivity deviation
  vector<Type> B_t = exp( log_B_t );
  vector<Type> Depletion_t = B_t / k;
  
  // Objective function
  vector<Type> nll_comp(5);
  nll_comp.setZero();

  recDev(0) = 0.;
  // Reconstruct time series
  for( int t=1; t<n_years; t++){
    Bpred_t[t]   = B_t[t-1]  + ((shape+1)/shape)*r*B_t[t-1]*(1-pow(sqrt(square((B_t[t-1]/k))),shape)) - c_t(t-1);
    recDev(t) = log(Bpred_t(t)/B_t(t));
    nll_comp(0) -= dnorm( log(B_t(t)), log(Bpred_t(t)), exp(log_sigmaP), true );
  }

  //calculate analytical q
  Type tmpSum=0.;
  Type tmpSumI=0.;
  Type tmpSumB=0.;
  Type nI=0; 
  for( int t=0; t<n_years; t++){
    if(!CppAD::isnan(log(I_t(t)))) nI += 1;
    if(!CppAD::isnan(log(I_t(t)))) tmpSum  += log(I_t(t) / B_t(t));
    if(!CppAD::isnan(log(I_t(t)))) tmpSumI += I_t(t);
    if(!CppAD::isnan(log(I_t(t)))) tmpSumB += B_t(t);
  }  
  Type q     = exp((1/nI) * tmpSum);
  

  //  CPUE likelihood
  for( int t=0; t<n_years; t++){
    if(!CppAD::isnan(log(I_t(t)))) nll_comp(1) -= dnorm( log(I_t(t)), log(q*B_t(t)), exp(log_sigmaI), true );   
  }                                                                                 
  
  // penalties for parameter priors
  nll_comp(2) -= dnorm( log_MSY, log(priorMode(1)), priorLogCV(1), true );
  nll_comp(2) -= dnorm( shape, priorMode(3), priorLogCV(3), true );
  nll_comp(2) -= dnorm( log_sigmaI, log(priorMode(4)), priorLogCV(4), true );
  nll_comp(2) -= dnorm( log_sigmaP, log(priorMode(5)), priorLogCV(5), true );

  // Penalty on starting biomass
  nll_comp(3) -= dnorm( log(Depletion_t(0)), log(priorMode(0)), priorLogCV(0), true );

  //try a flattened bottom U-shaped penalty function (platykurtic maybe)
  //value of boundWt likelihood units on bounds, rapidly increases    
  Type boundWt  = 1.;
  Type m        = (2*boundWt)/(log_kBounds(1) - log_kBounds(0));
  Type b        = boundWt - m*log_kBounds(1);
  Type boundPen = pow(m*log_k + b, 6);  
  nll_comp(4)  += boundPen;
  
  // Total likelihood
  Type nll = nll_comp.sum();

  // Reporting
  REPORT( log_B_t );
  REPORT( n_years );
  REPORT( B_t );
  REPORT( Bpred_t ); //not really useful - should convert to dev 
  REPORT( nll_comp );
  REPORT( recDev );
  REPORT( log_MSY );
  REPORT( k );
  REPORT( r );
  REPORT( shape );
  REPORT( q );
  REPORT( Depletion_t );
  REPORT( log_sigmaI );
  REPORT( log_sigmaP );

  ADREPORT( log_B_t );
  ADREPORT( B_t );
  ADREPORT( Depletion_t );
  ADREPORT( recDev );
  
  return nll;
}
