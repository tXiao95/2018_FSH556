#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // 
  DATA_VECTOR( y_t );
  DATA_IVECTOR( t );
  DATA_INTEGER( n );
  DATA_INTEGER( n_t );
  
  // prob of 0
  PARAMETER( theta );
  PARAMETER( log_sigma );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER( log_x0 );
  PARAMETER_VECTOR( log_x_t );
  
  Type p0 = exp(theta) / (1 + exp(theta));
  Type jnll=0;
  
  //X data latent
  jnll -= dnorm( log_x_t(0), alpha + rho*log_x0, exp(log_sigma), true );
  for( int j=1; j<n_t; j++ ){
    jnll -= dnorm( log_x_t(j),  alpha + rho*log_x_t(j-1), exp(log_sigma), true );
  }
  
  // Y data observed
  for( int i=0; i<n; i++ ){
    if( y_t(i)==0 ){
      jnll -= log(p0);
    } else{
      jnll -= log(1 - p0) + dnorm( log(y_t(i)), log_x_t(t(i)), exp(log_sigma), true );
    }
  }
  
  Type sigma = exp(log_sigma);
  
  REPORT( sigma );
  REPORT( p0 );
  REPORT( alpha );
  REPORT( rho );
  REPORT( log_x0 );
  REPORT( log_x_t );
  
  ADREPORT( sigma );
  ADREPORT( p0 );
  ADREPORT( alpha );
  ADREPORT( rho );
  ADREPORT( log_x0 );
  ADREPORT( log_x_t );
  
  return jnll;
}

