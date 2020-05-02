library(TMB)
library(data.table)
library(ggplot2)
library(lme4)
library(broom)

# Empirical Bayes random effect estimation

set.seed(1)
# Set up data -------------------------------------------------------------
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

# GLM: no variability -----------------------------------------------------
fit1 <- glm(y_i ~ 1, data=df, family = "poisson")

# GLM: Among-site variability ---------------------------------------------
fit2 <- glmer(y_i ~ (1 | group), data=df, family="poisson")

# TMB Models: overdispersion and among-site -----------------------------------------------------------
file <- "GLMM"
dyn.load( dynlib(file) )
nSim <- 100

# Indicate random effects (overdispersion and mean shift)
Random <- c("delta_i", "eps_s")

model <- 1:4

est <- CJ(model=1:4, nSim=1:100, mu=0, se=0)

for(i in 1:nSim){
  # Build inputs: Need size of arrays, and the design matrix itself
  df <- sim_data(1000, 10, 2, 1, sqrt(.5))
  Data <- list("n"=nrow(df), 
               "n_s"=length(unique(df$group)), 
               "y_i"=df$y_i, 
               "s_i"=df$group-1)
  
  #' Random: 1000 overdispersion terms, 10 group terms
  #' Fixed: 1 sigma_s, 1 sigma_y
  Parameters <- list( "log_sigma_s"=0, 
                      "log_sigma_y"=0, 
                      "mu"=0,
                      "delta_i"=rep(0, Data$n), 
                      "eps_s"=rep(0, Data$n_s))
  
  if(i %% 10 == 0) print(paste0("Sim ", i))
  for (m in model){
    Map <- list()
    #' Turn off overdispersion
    if(m %in% c(1,2)){
      Map[["log_sigma_y"]] <- factor(NA)
      Map[["delta_i"]] <- factor(rep(NA, Data$n))
    }
    
    #' Turn off site specific effect
    if(m %in% c(1,3)){
      Map[["log_sigma_s"]] <- factor(NA)
      Map[["eps_s"]] <- factor(rep(NA, Data$n_s))
    }
    Random_temp <- setdiff(Random, names(Map))
    Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random_temp, map=Map)
    # Turn off pritning of mgc (max gradient component)
    Obj$env$beSilent()
    Opt <- TMBhelper::fit_tmb(Obj)
    
    #' Matrix of the estimate and SE of each parameter estimate (random and fixed effects). 
    #' rownames are the parameters
    sd_matrix <- summary(sdreport(Obj))
    
    # Save estimate and SE
    est[model == m & nSim == i, `:=`(mu = sd_matrix['mu',1], 
                                     se = sd_matrix['mu', 2])]
  }
}

# Create mu table ---------------------------------------------------------
est[, model := ifelse(model == 1, "Global", 
                      ifelse(model == 2, "SiteEffect", 
                             ifelse(model == 3, "Overdispersion", "Full")))]

est[, `:=`(lower = mu - 1.96*se, upper = mu + 1.96*se)]
est[, containMu := 2 > lower & 2 < upper]

est[, .(mu_hat = mean(mu), coverageProb = mean(containMu)), model]

hist <- ggplot(est, aes(mu)) + geom_density(aes(col=model, fill=model), alpha=0.5) + 
  geom_vline(xintercept=2, linetype="dashed", col="red") + 
  theme_bw() + ggtitle("Density of mu-hat")

#' Density of mu estimates by simulation
pdf(file="hist.pdf", width=11, height=8.5)
print(hist)
dev.off()

#' CI coverage of truth
coverage <- ggplot(est, aes(nSim, mu)) + 
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper, col=containMu)) + 
  geom_hline(yintercept=2) + 
  facet_wrap(~model) + 
  theme_bw() + 
  ggtitle("Coverage")

pdf(file="coverage.pdf", width=11, height=8.5)
print(coverage)
dev.off()