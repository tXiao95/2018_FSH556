

library(TMB)
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 5 -- 1D spatial models/Lab/exercise/" )

###################
# Unequal distance autoregressive
###################

x = 1:100
Rho = 0.8
Sigma2 = (0.5) ^ 2
n_rep = 3
beta0 = 3

# Simulate locations
loc_s = rep(NA, length(x))
loc_s[1] = 0
for(s in 2:length(x)) loc_s[s] = loc_s[s-1] + rlnorm(1, meanlog=0, sdlog=1)

# Simulate spatial process
epsilon_s = rep(NA, length(x))
epsilon_s[1] = rnorm(1, mean=0, sd=sqrt(Sigma2))
for(s in 2:length(x)) epsilon_s[s] = Rho^abs(loc_s[s]-loc_s[s-1]) * epsilon_s[s-1] + rnorm(1, mean=0, sd=sqrt(Sigma2*(1-Rho^(2*abs(loc_s[s]-loc_s[s-1])))) )

# SImulate counts
c_si = matrix( nrow=length(x), ncol=n_rep)
for(s in 1:nrow(c_si)){
for(i in 1:ncol(c_si)){
  c_si[s,i] = rpois(1, exp(beta0 + epsilon_s[s]) )
}}

####################
# Fit model
####################

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_s"=rnorm(length(x)) )
compile( "Lab_5_exercise.cpp" )
dyn.load( dynlib("Lab_5_exercise") )

######## Version 2 -- Covariance and built-in function
# Build object
Data = list("Options_vec"=c(2), "c_si"=c_si, "loc_s"=loc_s )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="Lab_5_exercise" )
# Optimize
Opt2 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par2 = Opt2$par
h2 = Obj$env$spHess(random=TRUE)

######## Version 100 -- New Gaussian correlation function
# Build object
Data = list("Options_vec"=c(2), "c_si"=c_si, "loc_s"=loc_s )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="Lab_5_exercise" )
# Optimize
Opt100 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par100 = Opt100$par
h100 = Obj$env$spHess(random=TRUE)

# Compare timing
Opt2$run_time
Opt100$run_time

# Compare separability
library(INLA)
image(h2, main="Version 2"); dev.new()
image(h100, main="New option")

