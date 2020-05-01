#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Data integer array 
  DATA_INTEGER( n );
  DATA_INTEGER( n_s );
  // Can't be DATA_IARRAY if to be used in JNLL calculation
  DATA_VECTOR( y_i );
  DATA_IVECTOR( s_i );
  
  //Parameter scalars (fixed)
  PARAMETER( log_sigma_s );
  PARAMETER( log_sigma_y );
  PARAMETER( mu );
  
  //Parameter vectors (random)
  PARAMETER_VECTOR( delta_i );
  PARAMETER_VECTOR( eps_s );
  
  // variable to hold JNLL
  Type jnll = 0;
  
  //JNLL is composed of three parts
  
  // Y given delta_i and eps_s
  // Create temp variable to hold expected lambda
  vector<Type> ybar_i( n );
  for( int i=0; i<n; i++){
    ybar_i(i) = exp( mu + delta_i(i) + eps_s(s_i(i)) );
    jnll -= dpois( y_i(i), ybar_i(i), true );
  }
  
  // delta_i given sigma y (overdispersion)
  for( int i=0; i<n; i++){
    jnll -= dnorm( delta_i(i), Type(0.0), exp(log_sigma_y), true );
  }
  
  // eps_s given sigma_s (if among-site variability included)
  for( int s=0; s<n_s; s++){
    jnll -= dnorm( eps_s(s), Type(0.0), exp(log_sigma_s), true );
  }
  
  //Reporting
  Type sigma_y = exp(log_sigma_y);
  Type sigma_s = exp(log_sigma_s);
  
  REPORT( sigma_y );
  REPORT( sigma_s );
  REPORT( mu );
  
  ADREPORT( sigma_y );
  ADREPORT( sigma_s );
  ADREPORT( mu );
  
  return jnll;
}