library(TMB)
library(data.table)
library(ggplot2)

# Data sim ----------------------------------------------------------------

set.seed(10)
sim_data <- function(N, m, alpha, beta){
  #' N -- total number of observations
  #' m -- total number of groups
  #' alpha, beta: hyperparameters
   
  s_i <- sample(1:m, N, TRUE)
  p_s <- rbeta(m, alpha, beta)
  n_i <- sample(100:150, size=N, replace=TRUE)
  x_i <- sapply(1:N, function(i) rbinom(1, n_i[i], p_s[s_i[i]]))
  df  <- data.table(x_i=x_i, n_i = n_i, s_i=s_i, p = p_s[s_i])
  df
}

x <- seq(0,1,.01)
plot(x, dbeta(x, 50,50))

df <- sim_data(100, 5, .5, .5)
df[, p_hat := x_i / n_i]
df
plot(df$p, df$p_hat)

# TMB ---------------------------------------------------------------------
#' Fixed: alpha, beta
#' Random: p_s

Data <- list( "N"=nrow(df),
              "m"=length(unique(df$p)),
              "s_i"=df$s_i-1,
              "n_i"=df$n_i,
              "x_i"=df$x_i )

Parameters <- list( "log_alpha"=0, 
                    "log_beta"=0,
                    "p_s"=rep(0.5, Data$m))

Random <- c("p_s")

compile( "betaBinomial.cpp" )
dyn.load( dynlib("betaBinomial"))
Obj <- MakeADFun(Data, Parameters, random=Random)

Obj$env$beSilent()

Opt <- TMBhelper::fit_tmb(Obj)

# simulation --------------------------------------------------------------
nObs <- 100
nGroup <- 5
nSim <- 100

# Breaks when prior for p is not unimodal 
a <- c(.5, 5, 1, 2, 2)
b <- c(.5, 1, 3, 2, 5)

a <- c(2)
b <- c(2)
nPar <- length(a)

est <- data.table(log_alpha=log(a), log_beta=log(b), mean=0, se=0)
#est <- data.table(param=c("log_alpha", "log_beta"), mean=0, se=0)
#est <- est[rep(1:.N, nPar)]
#est[param == "log_alpha", truth := log(a)]
#est[param == "log_beta", truth := log(b)]
est <- est[rep(1:.N, nSim)][, sim := 1:.N, by=.(log_alpha, log_beta)]
est[, param:=.GRP, by=.(log_alpha, log_beta)]

for(i in 1:nSim){
  message(i)
  for(j in 1:nPar){
    df <- sim_data(nObs, nGroup, a[j], b[j])
    Data <- list( "N"=nrow(df),
                  "m"=length(unique(df$p)),
                  "s_i"=df$s_i-1,
                  "n_i"=df$n_i,
                  "x_i"=df$x_i )
    
    Parameters <- list( "log_alpha"=0, 
                        "log_beta"=0,
                        "p_s"=rep(0.5, Data$m))
    Random <- c("p_s")
    
    Obj <- MakeADFun(Data, Parameters, random=Random)
    Obj$env$beSilent()
    Opt <- TMBhelper::fit_tmb(Obj)
    
    est_matrix <- summary(Opt$SD)
    
    b_est <- est_matrix['log_beta', 1]
    b_se <- est_matrix['log_beta', 2]
    a_est <- est_matrix['log_alpha', 1]
    a_se <- est_matrix['log_beta', 2]
    
    est[sim==i & log_alpha==log(a[j]) & log_beta==log(b[j]), `:=`(log_alpha_est=a_est, log_alpha_se=a_se,
                                                log_beta_est=b_est, log_beta_se=b_se)]
  }
}

est[, lower_log_alpha := log_alpha_est - 1.96*log_alpha_se]
est[, upper_log_alpha := log_alpha_est + 1.96*log_alpha_se]
est[, contain_alpha := log_alpha > lower_log_alpha & log_alpha < upper_log_alpha]

est[, lower_log_beta := log_beta_est - 1.96*log_beta_se]
est[, upper_log_beta := log_beta_est + 1.96*log_beta_se]
est[, contain_beta := log_beta > lower_log_beta & log_beta < upper_log_beta]

est[, mean(contain_alpha), param]
est[, mean(contain_beta), param]

ggplot(est, aes(sim, log_alpha_est)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_log_alpha, ymax=upper_log_alpha, col=contain_alpha)) + 
  facet_wrap(~param) + 
  geom_hline(aes(yintercept=log_alpha))

ggplot(est, aes(sim, log_beta_est)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_log_beta, ymax=upper_log_beta, col=contain_beta)) + 
  facet_wrap(~param) + 
  geom_hline(aes(yintercept=log_beta))
