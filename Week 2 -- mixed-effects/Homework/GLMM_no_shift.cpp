#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Data integer array 
  DATA_INTEGER( n );
  // Can't be DATA_IARRAY if to be used in JNLL calculation
  DATA_VECTOR( y_i );
  
  //Parameter scalars (fixed)
  PARAMETER( log_sigma_y );
  PARAMETER( mu );
  
  //Parameter vectors (random)
  PARAMETER_VECTOR( delta_i );
  
  // variable to hold JNLL
  Type jnll = 0;
  
  //JNLL is composed of three parts
  
  // Y given delta_i and eps_s
  // Create temp variable to hold expected lambda
  vector<Type> ybar_i( n );
  for( int i=0; i<n; i++){
    ybar_i(i) = exp( mu + delta_i(i) );
    jnll -= dpois( y_i(i), ybar_i(i), true );
  }
  
  // delta_i given sigma y (overdispersion)
  for( int i=0; i<n; i++){
    jnll -= dnorm( delta_i(i), Type(0.0), exp(log_sigma_y), true );
  }
  
  //Reporting
  Type sigma_y = exp(log_sigma_y);
  
  REPORT( sigma_y );
  REPORT( mu );
  
  ADREPORT( sigma_y );
  ADREPORT( mu );
  
  return jnll;
}