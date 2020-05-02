#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER( N );
  DATA_INTEGER( m );
  // This is for storing indices so okay that it doesn't take 
  // DATA_VECTOR of type Type.
  DATA_IVECTOR( s_i );
  // Even though these are integer data types, make sure
  // the C++ functions they are passed to take the right type. 
  // dbeta and dbinom take 'Type' types
  DATA_VECTOR( n_i );
  DATA_VECTOR( x_i );
  
  PARAMETER( log_alpha );
  PARAMETER( log_beta );
  PARAMETER_VECTOR( p_s );
  
  Type jnll=0;
  
  // Data likelihood to p
  for( int i=0; i<N; i++){
    jnll -= dbinom( x_i(i), n_i(i), p_s(s_i(i)), true );
  }
  
  // p likelihood to alpha beta
  for( int j=0; j<m; j++){
    jnll -= dbeta( p_s(j), exp(log_alpha), exp(log_beta), true );
  }
  
  Type alpha = exp(log_alpha);
  Type beta = exp(log_beta);
  
  REPORT( log_alpha );
  REPORT( log_beta );
  REPORT( alpha );
  REPORT( beta );
  REPORT( p_s );
  
  ADREPORT( log_alpha );
  ADREPORT( log_beta );
  ADREPORT( p_s );
  
  return jnll;
}