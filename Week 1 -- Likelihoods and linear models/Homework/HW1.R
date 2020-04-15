library(data.table)
library(ggplot2)
library(here)
library(statmod)
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
compile( "deltaModels.cpp" )
catch <- EBS$catch
lat <- EBS$catch
X = cbind( "Intercept"=rep(1,length(catch)))
k_i = rep(1, length(catch))

# Step 2 -- build inputs and object
dyn.load( dynlib("deltaModels") )
Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))
Data = list( "y_i"=catch, "X_ij"=X, "Options_vec"=c(0), "k_i"=k_i, k=c(2))
Obj = MakeADFun( data=Data, parameters=Params, DLL="deltaModels")

# Step 3 -- test and optimize
initial_loss <- Obj$fn( Obj$par )
initial_grad <- Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), 
                              "Est"=Opt$par, 
                              "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors

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
sprintf("SD of log catch from model %f", (Opt$par[3]))

# Cross validation --------------------------------------------------------
nFolds <- 10
EBS[, k := sample(1:nFolds, size=.N, replace = TRUE)]
TMBfile <- "deltaModels"

compile( sprintf("%s.cpp", TMBfile) )
dyn.load( dynlib(TMBfile) )
catch <- EBS$catch
lat   <- EBS$catch

Data = list( "y_i"=catch, "X_ij"=X, "k_i"=EBS$k )
X = cbind( "Intercept"=rep(1,length(catch)))
Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))

loss <- data.table::CJ(k = 1:nFolds,
                   opt = 0:2,
                   train_loss = 0, 
                   test_loss = 0 )

# For each likelihood, do 10-fold CV
for (option in 0:2){
  for(fold in 1:nFolds){
    sprintf("Option %d, fold %d", option, fold)
    Data$Options_vec <- option
    Data$k <- fold
    
    nTest <- nrow(EBS[k == fold])
    nTrain <- nrow(EBS[k != fold])
    
    Obj = MakeADFun( data=Data, parameters=Params, DLL="deltaModels")
    # initial_loss <- Obj$fn( Obj$par )
    # initial_grad <- Obj$gr( Obj$par )
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
    Opt$diagnostics = data.frame( "name"=names(Obj$par), 
                                  "Est"=Opt$par, 
                                  "final_gradient"=as.vector(Obj$gr(Opt$par)))
    # Opt$par # estimated parameters
    SD = sdreport( Obj ) # standard errors
    loss[k == fold & opt == option, `:=`(train_loss = Obj$report()$train_jnll / nTrain, 
                                   test_loss = Obj$report()$test_jnll / nTest)]
  }
}
# negative log likelihood
Obj$fn( Opt$par )
# Gradient at the last iteration
Obj$gr( Opt$par )
# Log preditive score per datum


loss[, mean(test_loss), by = opt]
# Seems like the delta-lognormal model has the lowest out of sample error

# Simulation --------------------------------------------------------------

#' The following experiment is as follows. We generate some simulated data based on 
#' three different probability distributions. They are
#' 
#' 1) Lognormal
#' 2) Gamma
#' 3) Inverse gaussian
#' 
#' Then we estimate the parameters of interest using three different estimators. These are the
#' three different models we estimated above. We perform a 3 x 3 analysis to estimate
#' paramters under each unique scenario. 
#' 
#' 

# Simulate data. For the distributions we simulate data from, we fix the parameters
#' at what the models estimated from the entire dataset. 
EBS$k <- -999
X = cbind( "Intercept"=rep(1,length(EBS$catch)))
# I have data arguments for K_i (the fold the row is in) and k (the fold I'm testing on). 
#' Since I just want to train on everything, then just pick random values so that don't overlap
Data = list( "y_i"=EBS$catch, "X_ij"=X, "k_i"=EBS$k, "k" = 999)
Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))

TMBfile <- "deltaModels"
dyn.load( dynlib(TMBfile) )

n <- nrow(EBS)
nsim <- 100

data <- CJ(sim=1:nsim, simM = 0:2, modelM = 0:2)

# Hopefully these work like numpy?
truth <- array(0, c(3, 3))
simData <- array(0, c(3, 100, 3, 3))
# 3 sim models, 100 replications, 3 estimation methods, 3 parameters each

for(simModel in 1:3){
  # Estimate a model on the whole dataset using 
  # 1: lognormal, 2: gamma, 3: inverse gaussian
  Data = list( "y_i"=EBS$catch, "X_ij"=X, "k_i"=EBS$k, "k" = 999, "Options_vec"=c(simModel-1))
  Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))
  ObjTruth = MakeADFun( data=Data, parameters=Params, DLL="deltaModels")
  OptTruth = nlminb( start=ObjTruth$par, objective=ObjTruth$fn, gradient=ObjTruth$gr )
  truth[simModel,] <- OptTruth$par

  #' nsim iterations to get a distribution around parameters
  for(i in 1:nsim){
    p0 <- exp(OptTruth$par[2])
    if(simModel == 1){
      # Lognormal
      logmean <- OptTruth$par[1]
      logsd <- OptTruth$par[3]
      simY <- rbinom(n, 1, 1-p0) * rlnorm(n, logmean, logsd)
    } else if(simModel == 2){
      # Gamma
      xb <- OptTruth$par[1]
      k <- 1 / OptTruth$par[3]^2
      theta <- exp(xb) / k
      simY <- rbinom(n, 1, 1-p0) * rgamma(n, shape=k, scale=theta)
    } else if(simModel == 3){
      # Inverse Gaussian
      mean <- OptTruth$par[1]
      shape <- OptTruth$par[3]
      simY <- rbinom(n, 1, 1-p0) * rinvgauss(n, exp(mean), shape)
    }
    for(modelNum in 1:3){
      #' Use each estimation strategy on each data generating process
      Data = list( "y_i"=simY, "X_ij"=X, "k_i"=EBS$k, "k" = 999, "Options_vec" = modelNum-1)
      Params = list("b_j"=rep(1,ncol(X)), "theta_z"=c(1,1))
      ObjSim = MakeADFun( data=Data, parameters=Params, DLL="deltaModels")
      OptSim = nlminb( start=ObjSim$par, objective=ObjSim$fn, gradient=ObjSim$gr )
      
      #' Store parameter estimates in array
      simData[simModel, i, modelNum, ] <- OptSim$par
    }
  }
}

# Collapse to mean information
simMeans <- apply(simData, c(1, 3, 4), mean)

# Plotting simulation results ---------------------------------------------

# Lognormal Simulation results
simMeans[1,1,]
truth[1,]

# Gamma simulation results
simMeans[2,2,]
truth[2,]

# Invgauss simulation resutls
simMeans[3,3,]
truth[3,]

#' Plots
#' For each simulation model
#' 
#' Facet by (Estimat x parameter). Overlaid with the truth for each one

key <- CJ(modelNum=1:3, paramNum=1:3)
key[, modelName := ifelse(modelNum == 1, "LN", 
                          ifelse(modelNum == 2, "Gamma", "IG"))]
key$param <- c("E(log(Y))", "p0", "SE(log(Y))", 
               "log(E(Y))", "p0", "CV", 
               "log(E(Y))", "p0", "shape")

pdf(file="simulations.pdf", width=11, height=8.5)
for(simNum in 1:3){
  simModelName <- key[modelNum == simNum, unique(modelName)]
  par(mfrow=c(3, 3))
  for(estNum in 1:3){
    for(parameter in 1:3){
      keysub <- key[modelNum == estNum & paramNum == parameter]
      data <- simData[simNum,,estNum,parameter]
      hist(data, main=sprintf("%s Est for %s", keysub[,modelName], keysub[,param]))
      abline(v=truth[estNum, parameter], col="red", lwd=3, lty=2)
    }
  }
  mtext(sprintf("%s Simulated", simModelName), line = -1.4, outer = TRUE)
}
dev.off()
