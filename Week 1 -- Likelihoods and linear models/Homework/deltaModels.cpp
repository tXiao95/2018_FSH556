
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  // I'm not sure why originally this code had subtracted the log(x). It's not like we worry
  // about underflow or overflow with the likelihood itself, but rather the exp term within dnorm.
  //Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  Type logres = dnorm( log(x), meanlog, sdlog, true);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Integers give settings
  // Options_vec: Whether lognormal, gamma, invgauss
  // k: Which fold to test on
  DATA_IVECTOR( Options_vec );
  DATA_IVECTOR( k );
  
  // Fold indicator, response and design matrix
  DATA_VECTOR( k_i );
  DATA_VECTOR( y_i );
  DATA_MATRIX( X_ij );
  
  // Parameters
  PARAMETER_VECTOR( b_j );
  //  theta[0] = P(Y=0), theta[1] = sd of Y, Y>0
  PARAMETER_VECTOR( theta_z );
  
  // Logistic model for Y=0, Lognormal model for Y>0, SD of log(Y)
  Type zero_prob = 1 / (1 + exp(-theta_z(0)));
  Type logsd = theta_z(1);
  
  // Initialize cost: joint negative log likelihood
  Type jnll = 0;
  Type train_jnll = 0;
  Type test_jnll = 0;
  
  // Number of observations, N
  int n_data = y_i.size();
  
  // Linear predictor - X*B, (n x p) * (p x 1)
  vector<Type> linpred_i( n_data );
  linpred_i = X_ij * b_j;
  
  // Loop through each data point to calculate additive contribution to JNLL
  for( int i=0; i<n_data; i++){
    if(y_i(i)==0){
      jnll = log( zero_prob );
    } else{
      if( Options_vec(0)==0 ) jnll = log( 1-zero_prob ) + dlognorm( y_i(i), linpred_i(i), logsd, true );
      if( Options_vec(0)==1 ) jnll = log( 1-zero_prob ) + dgamma( y_i(i), 1/pow(logsd,2), exp(linpred_i(i))*pow(logsd,2), true );
      if( Options_vec(0)==2 ) jnll = log( 1-zero_prob ) + dinvgauss( y_i(i), exp(linpred_i(i)), logsd, true );
    }
    if(k_i(i)==k(0)){
      test_jnll -= jnll;
    } else{
      train_jnll -= jnll;
    }
  }
  
  // Reporting
  REPORT( zero_prob );
  REPORT( logsd );
  REPORT( linpred_i );
  ADREPORT( zero_prob );
  ADREPORT( logsd );
  REPORT( test_jnll );
  REPORT( train_jnll );
  
  return train_jnll;
}
