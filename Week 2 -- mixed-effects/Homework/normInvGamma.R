library(TMB)
library(data.table)
library(MCMCpack)
library(ggplot2)

# Data Generation ---------------------------------------------------------

#' Empirical Bayes estimation of Normal likelihood, normal prior on mu, log-normal prior on sigma

sim_data <- function(n, J, mu0, sigma_mu, l_mu_sigma, l_sigma_sigma){
  #' sigma_s   ~ LogN(l_mu_sigma, l_sigma_sigma)
  #' mu_s      ~ N(mu_0, sigma_mu)
  #' Y         ~ N(mu_s, sigma_s)
  #' 
  #' Fixed: mu_0, sigma_mu, l_mu_sigma, l_sigma_sigma
  #' Random: sigma_s, mu_s
  
  #' Group labels of the data  
  s_i <- sample(1:J, n, replace=TRUE)
  
  #' site specific means
  mu_s <- rnorm(J, mu0, sigma_mu)
  
  #' site specific variance
  log_sigma_s <- rnorm(J, l_mu_sigma, l_sigma_sigma)
  
  #' Observed data
  y_i <- sapply(s_i, function(j) rnorm(1, mu_s[j], exp(log_sigma_s[j])))
  
  df <- data.table(y_i=y_i, s_i=s_i)
  return(df)
}

df <- sim_data(1000, 10, 10, 3, 1, .5)

ggplot(df, aes(as.factor(s_i), y_i)) + geom_boxplot()

# TMB ---------------------------------------------------------------------

#' TMB only works well with normal priors...
#' Wanted to do InvGamma on variance, but maybe will try LogNormal

Data <- list( "n" = nrow(df), 
              "n_s" = length(unique(df$s_i)), 
              "y_i" = df$y_i, 
              "s_i" = df$s_i-1 )

#' l_mu_sigma: mean of log(sigma_s)
#' l_sigma_sigma: sigma of log(sigma_s)
Parameters <- list( "mu_0" = 0,
                    "log_sigma_mu" = 0,
                    "l_mu_sigma" = 0, 
                    "l_sigma_sigma" = 1, 
                    "mu_s" = rep(0, Data$n_s), 
                    "log_sigma_s" = rep(0, Data$n_s))

Random <- c("mu_s", "log_sigma_s")

compile("normInvGamma.cpp")
dyn.load( dynlib("normInvGamma"))

Obj <- MakeADFun(data = Data, parameters = Parameters, random = Random)
Obj$env$beSilent()
Opt <- TMBhelper::fit_tmb(Obj)

sd_table <- summary(sdreport(Obj))[c('mu_0', 'sigma_mu', 'l_mu_sigma', 'l_sigma_sigma'), ]
sd_table
