library(TMB)
library(data.table)
library(ggplot2)
library(lme4)
library(broom)

set.seed(1)
# Set up data -------------------------------------------------------------
n <- 1000
n_s <- 10
mu <- 2

sim_data <- function(n, n_s, mu, sigma_s, sigma_y){
  # Group level effect
  eps_s <- rnorm(n_s, 0, sigma_s)
  group <- sample(1:n_s, n, replace=TRUE)
  mean_s <- eps_s[group]
  
  # Overdispersion effect
  delta_i <- rnorm(n, 0, sigma_y)
  
  exp_count <- exp(mu + mean_s + delta_i)
  y_i <- sapply(exp_count, function(lambda) rpois(1, lambda))
  data.table(y_i=y_i, group=group)
}

df <- sim_data(1000, 10, 2, 1, sqrt(.5))

ggplot(df, aes(as.factor(group), y_i)) + geom_boxplot()

# GLM: no variability -----------------------------------------------------
fit1 <- glm(y_i ~ 1, data=df, family = "poisson")

# GLM: Among-site variability ---------------------------------------------
fit2 <- glmer(y_i ~ (1 | as.factor(group)), data=df, family="poisson")

# TMB: just overdispersion ------------------------------------------------
file_no_shift <- "GLMM_no_shift"
compile(paste0(file_no_shift, ".cpp"))

Data <- list("n"=nrow(df), 
             "y_i"=df$y_i)

Parameters <- list( "log_sigma_y"=0,
                    "mu"=0,
                    "delta_i"=rep(0, Data$n))
Random <- c("delta_i")

dyn.load( dynlib(file_no_shift) )
Obj_no_shift <- MakeADFun(Data, Parameters, random = Random)
Opt_no_shift <- TMBhelper::fit_tmb(Obj_no_shift)

# Fit TMB Model: overdispersion and among-site -----------------------------------------------------------
file <- "GLMM"
compile(paste0(file, ".cpp"))

# Build inputs: Need size of arrays, and the design matrix itself
Data <- list("n"=nrow(df), 
             "n_s"=length(unique(df$group)), 
             "y_i"=df$y_i, 
             "s_i"=df$group-1)

#' Parameters to estimate: 
#' Random: 1000 overdispersion terms, 10 group terms
#' Fixed: 1 sigma_s, 1 sigma_y
Parameters = list( "log_sigma_s"=0, 
                   "log_sigma_y"=0, 
                   "mu"=0,
                   "delta_i"=rep(0, Data$n), 
                   "eps_s"=rep(0, Data$n_s))

# Indicate random effects (overdispersion and mean shift)
Random = c("delta_i", "eps_s")

# Load .cpp template
dyn.load( dynlib(file) )

# Create TMB object
Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random)

# Minimize
nlOpt <- nlminb(Obj$par, Obj$fn, Obj$gr)
Opt <- TMBhelper::fit_tmb(Obj)


# Create mu table ---------------------------------------------------------

mu_table <- data.table(model=c("GLM", "GLMER", "TMB_no_shift", "TMB_full"))

mu_table[, mu := c(fit1$coefficients, 
                   tidy(fit2)$estimate[1], 
                   Opt_no_shift$SD$value[2], 
                   Opt$SD$value[3])]

mu_table[, sd := c(tidy(fit1)$std.error, 
                   tidy(fit2)$std.error[1], 
                   Opt_no_shift$SD$sd[2], 
                   Opt$SD$sd[3])]

mtidy(fit1)
tidy(fit2)

Obj$report()

Opt
Opt_no_shift
