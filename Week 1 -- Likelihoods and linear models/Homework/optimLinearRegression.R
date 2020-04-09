
###### NORMAL #####################################
# Objective function for iid normal data
obj <- function(params, x, y){
  beta0 <- params[1]
  beta1 <- params[2]
  sigma <- params[3]
  # Negative log likelihood to minimize
  -sum(log(dnorm(y, beta0 + beta1*x, sigma)))
}

# b0 = 2, b1 = 3, sigma = 1
x <- rnorm(100)
y <- 3*x + 2 + rnorm(100, 0, 1)

optimizer <- optim(par=c(0, 0, 1), fn=obj, x=x, y=y)
nl_opt <- nlminb(start=c(0, 0, 1), objective=obj, y=y, x=x)

optimizer
nl_opt

####### BERNOULLI ####################################
# Objective function for iid bernoulli data (logit link)
obj2 <- function(params, x, y){
  beta0 <- params[1]
  beta1 <- params[2]
  lin_pred <- beta1*x + beta0
  p <- exp(lin_pred) / (1 + exp(lin_pred))
  -sum(y*log(p) + (1-y)*log(1-p))
}

# b0 = 1, b1 = 1
x <- rnorm(100)
l <- x + 1; p <- exp(l) / (1 + exp(l))
y <- rbinom(100, 1, p)

optimizer <- optim(par=c(0, 0), fn=obj2, x=x, y=y)
nl_opt <- nlminb(start=c(0, 0), objective=obj2, y=y, x=x)
optimizer
nl_opt
