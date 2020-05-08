library(FishStatsUtils)
library(data.table)
library(TMB)
library(ggplot2)

# Setup -------------------------------------------------------------------

data("EBS_pollock_data")
df <- data.table(EBS_pollock_data)
df[, t := year - 1982]

Data <- list( y_t = df$catch, 
              t = df$t, 
              n = nrow(df), 
              n_t = length(unique(df$t)) )

Parameters <- list( theta = 0, 
                    log_sigma = 0,
                    alpha = 0, 
                    rho = 0, 
                    log_x0 = 0,
                    log_x_t = rep(1, Data$n_t) )

Random <- c("log_x_t")

compile( "state_space.cpp" )
dyn.load( dynlib("state_space") )

Obj <- MakeADFun(Data, Parameters, random=Random)

Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt <- TMBhelper::fit_tmb(Obj)

# Results -----------------------------------------------------------------

estnames <- rownames(summary(Opt$SD))
est <- data.table(summary(Opt$SD))
est[, param := estnames]
est <- est[1:38]

est_x <- est[6:nrow(est),]
est_x[, sigma := exp(est[param == "log_sigma", Estimate])]

est_x[, year:= as.factor(unique(df$year))]
est_x[, lower := Estimate - 1.96*sigma]
est_x[, upper := Estimate + 1.96*sigma]
est_x[, real_mean := exp(Estimate + sigma^2/2)]
est_x[, real_sigma := sqrt((exp(sigma^2)-1)*exp(2*Estimate+sigma^2))]

density <- CJ(year=as.factor(1982:2014), catch = seq(-10, 10, .1))
density <- merge(density, est_x[, .(year, Estimate, sigma)], by="year")
density[, y:=dnorm(catch, Estimate, sigma)]

boxplots <- ggplot() + 
  geom_boxplot(data=df[catch!=0], aes(as.factor(year), log(catch))) + 
  geom_point(data=est_x, aes(year, Estimate), col="red") + 
  geom_errorbar(data=est_x, aes(year, Estimate, ymin = lower, ymax = upper), col="red", alpha = 0.5, size=2)

density_plot <- ggplot(df, aes(log(catch))) + 
  geom_density(aes(col="data")) + 
  facet_wrap(~as.factor(year)) + 
  theme_bw() + 
  geom_line(data=density, aes(catch, y, col="fit"))

pdf(file="density_of_fit.pdf", width=11, height=8.5)
print(density_plot)
dev.off()

pdf(file="boxplot_of_fit.pdf", width=11, height=8.5)
print(boxplots)
dev.off()
