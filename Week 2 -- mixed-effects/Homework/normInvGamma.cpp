#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER( n );
  DATA_INTEGER( n_s );
  DATA_VECTOR( y_i );
  DATA_IVECTOR( s_i );
  
  //Fixed effects estimate by EB
  PARAMETER( mu_0 );
  PARAMETER( log_sigma_mu );
  PARAMETER( l_mu_sigma );
  PARAMETER( l_sigma_sigma );
  
  //Random effects
  PARAMETER_VECTOR( mu_s );
  PARAMETER_VECTOR( log_sigma_s );
  
  Type jnll=0;
  
  // y_i data likelihood
  for(int i=0; i<n; i++) {
    jnll -= dnorm( y_i(i), mu_s(s_i(i)), exp( log_sigma_s(s_i(i)) ), true );
  }
  
  // mu_s random likelihood
  for(int j=0; j<n_s; j++){
    jnll -= dnorm( mu_s(j), mu_0, exp(log_sigma_mu), true );
  }
  
  // log_sigma_s random likelihood
  for(int j=0; j<n_s; j++){
    jnll -= dnorm( log_sigma_s(j), l_mu_sigma, l_sigma_sigma, true );
  }
  
  Type sigma_mu = exp(log_sigma_mu);
  vector<Type> sigma_s = exp(log_sigma_s);
  
  ADREPORT( mu_0 );
  ADREPORT( sigma_mu );
  ADREPORT( l_mu_sigma );
  ADREPORT( l_sigma_sigma );
  ADREPORT( mu_s );
  ADREPORT( sigma_s );
  
  REPORT( mu_0 );
  REPORT( sigma_mu );
  REPORT( l_mu_sigma );
  REPORT( l_sigma_sigma );
  REPORT( mu_s );
  REPORT( sigma_s );
  
  return jnll;
}