library(data.table)
library(ggplot2)
library(here)
library(TMB)
library(VAST)

data( EBS_pollock_data, package="FishStatsUtils")
EBS <- data.table(EBS_pollock_data)
propzero <- EBS[, .N, .(catch>0, year)][, rate := N / sum(N)][catch == FALSE]
zeros <- EBS[catch == 0]
zeros_summary <- zeros[,.N,year]

# Exploratory -------------------------------------------------------------

summary(EBS)
#' My understanding is that the 'catch' is the number of fish caught per 
#' fishing trip. So each row of this dataset is what unit exactly? The numbers
#' are also not integers, so it's not number of fish. There is one observation
#' per (lat, long, year).

# 'catch' indicates catch rate of the pollock
ggplot(EBS, aes(catch)) + 
  geom_histogram(bins=100) + 
  ggtitle("Distribution catch rate")
#' 424 zeros in the data. Some said fishers
EBS[, .N, catch == 0]
#' A log transform on the non zero data show skewed left distribution of catch rate. 
#' One could maybe try a log normal, that would underestimate lower values and overestimate
#' extremely large values in this distribution. 
ggplot(EBS[catch != 0], aes(catch)) + geom_histogram() + scale_x_log10()

#' Plots number of locations that had 0 fish caught. That spike from 2007-2010 is interesting. 
ggplot(zeros_summary, aes(year, N)) + 
  geom_line() + 
  geom_point()

#' Doesn't add much information since the number of lat-longs in each year is roughly 350-376 
#' for every year
ggplot(propzero, aes(year, rate)) + 
  geom_line() + 
  geom_point()

#' Distribution of log catch rate. The zeros show a loss in fish catching starting 2007-2010 and
#' a bounceback in 2011. But in the non zero catch data, the decline in median catch starts in 2003
ggplot(EBS, aes(as.factor(year), catch)) + 
  geom_boxplot() + 
  scale_y_log10()

EBS[, logcatch := log(catch)]
logEBS <- EBS[catch != 0, .(med=median(logcatch), 
        lower = quantile(logcatch, .025), 
        upper = quantile(logcatch, .975)), year]

#' Can see the dip in catch rates starts as early as 2004. In 2004, median is 28.65
#' fish, then it's as low as 1.84 fish in 2009. Is it related to the Great Recession? Huge bounceback
#' though. It seems much more likely to do with the economy and less the natural world. Should see if
#' the water temperature covariates have anything to do with it. 
ggplot(logEBS, aes(year, med)) + 
  geom_line() + 
  geom_point() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.2)

# Model Fitting -----------------------------------------------------------
 # Compile code
compile( "deltaLogNormal.cpp" )
catch <- EBS$catch
lat <- EBS$catch
X = cbind( "Intercept"=rep(1,length(catch)))

# Step 2 -- build inputs and object
dyn.load( dynlib("deltaLogNormal") )
Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))
Data = list( "y_i"=catch, "X_ij"=X, "Options_vec"=c(0))
Obj = MakeADFun( data=Data, parameters=Params, DLL="deltaLogNormal")

# Step 3 -- test and optimize
initial_loss <- Obj$fn( Obj$par )
initial_grad <- Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), 
                              "Est"=Opt$par, 
                              "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors

### Check convergence
# Are the gradients close to zero
all(abs(Opt$diagnostics[,'final_gradient'])<0.001)
# How much did we reduce the gradient from the start by a factor of?
norm(initial_grad, "2") / norm(Opt$diagnostics$final_gradient, "2")
# Pretty large, so even though we don't have exactly 0 gradients, we've come
# Is the hessian positive definite?
(SD$pdHess==TRUE)

# Compare estimates to data for Lognormal -----------------------------------------------

prob0 <- length(catch[catch==0]) / length(catch)
sprintf("Probability of 0 from data: %f", prob0)
sprintf("Probability of 0 from model: %f", exp(Opt$par[2]))
logmean <- mean(log(catch[catch > 0]))
sprintf("Mean of log catch from data %f", logmean)
sprintf("Mean of log catch from model %f", (Opt$par[1]))
sprintf("Median of catch from data %f", exp(logmean))
sprintf("Exp(mu) of catch from model %f", exp(Opt$par[1]))
logsd <- sd(log(catch[catch > 0]))
sprintf("SD of log catch from data %f", logsd)
sprintf("SD of log catch from model %f", exp(Opt$par[3]))

