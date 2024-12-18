library(mvtnorm)
set.seed(571)
rm(list = ls())
alpha=0.05
N <- 1e2

# ~~~~~~~~~ #
# Problem 1 #
# ~~~~~~~~~ #
m=5; n=4; mu=0; rho=.9; beta0=1; beta1=2.5; sigma=1

simul_reg <- function(m, n, mu=0, rho, beta0=1, beta1=0,
                      sigma=1, alpha=0.05) {
  
  # build varcovar matrix
  # random intercept is N(0, theta)
  # corr(Y_ij, Y_ik) = theta/(theta * sigma)
  # where sigma is the variance of the truly stochastic error e~N(0,sigma^2)
  # we set corr(Y_ij, Y_ik) = to rho, then we can solve for theta
  # theta = (rho*sigma^2)/(rho-1)
  theta = -(rho*sigma^2)/(rho-1)
  varY <- theta + sigma^2
  covY <- rho*varY
  varcovar <- matrix(covY, nrow=n, ncol=n)
  diag(varcovar) <- varY
  
  # Create correlated e vector
  e <- mvtnorm::rmvnorm(n=m, mean=rep(mu,n), sigma=varcovar) |>
       matrix(ncol=1, nrow=m*n)
  id <- rep(1:m,n)
  
  # Create X vector where trt is randomly assigned per id
  trt <- sample(c(0,1),1)
  X <- rep(trt, m*n)
  if(trt==1) { 
    X[id %% 2 == 0] <- 0 
  }else { 
    X[id %% 2 == 0] <- 1 
  }
  
  # Create Y vector with true relationship
  Y <- beta0 + beta1*X + e 
  
  # Fit linear model
  fit <- lm(Y~X)
  beta_hat <- coef(fit)
  beta_hat_SE <- sqrt(diag(vcov(fit)))
  CI_beta_hat <- c(beta_hat - abs(qt(alpha/2, m*n-2))*beta_hat_SE, 
                   beta_hat + abs(qt(alpha/2, m*n-2))*beta_hat_SE)
  inside <- c(as.numeric(CI_beta_hat[[1]] <= beta0 & beta0 <= CI_beta_hat[[3]]),
              as.numeric(CI_beta_hat[[2]] <= beta1 & beta1 <= CI_beta_hat[[4]]))
  pval <- summary(fit)$coefficients[,4]

  #type I error: rejecting the null when the null is true
  if(beta1==0 | beta0==0) {
    type_1 <- c(as.numeric(pval[[1]]==1 & beta0==0),
                as.numeric(pval[[2]]==1 & beta1==0))
  } else{
    type_1 <- c(NA,NA)
  }
  
  o <- list("b0_inside" = inside[[1]], 
            "b0_type_1" = type_1[[1]], 
            "b0_hat"    = beta_hat[[1]],
            "b1_inside" = inside[[2]], 
            "b1_type_1" = type_1[[2]], 
            "b1_hat"    = beta_hat[[2]])
  return(o)
}


#################### create table of results
run_sims <- function(N, m, n, mu, rho, beta0, beta1, alpha=0.05, sigma=1) {
  coverage_b0 <- rep(0, N)   ;coverage_b1 <- rep(0, N)
  type_1_b0   <- rep(0,N)    ;type_1_b1   <- rep(0,N)
  b0_hat      <- rep(0,N)    ;b1_hat      <- rep(0,N)
  for(i in 1:N) {
    r <- simul_reg(m, n, mu, rho, beta0, beta1, alpha, sigma)
    coverage_b0[i] <- r[[1]] 
    type_1_b0[i] <- r[[2]] 
    b0_hat[i] <- r[[3]]
    coverage_b1[i] <- r[[4]] 
    type_1_b1[i] <- r[[5]] 
    b1_hat[i] <- r[[6]]    
  }
  o <- list("coverage_b0" = mean(coverage_b0) ,"coverage_b1" = mean(coverage_b1),
            "type_1_b0"   = mean(type_1_b0)   ,"type_1_b1"   = mean(type_1_b1),
            "b0_hat"      = beta0-mean(b0_hat),"b1_hat"      = beta1-mean(b1_hat))
  return(o)
}

run_sims(N=N, m=m, n=n, mu=0, rho=rho, beta1=0, beta0=0)

tab <- dplyr::tibble(beta1=numeric(), beta0=numeric(),
                     rho=numeric(), m=numeric(), n=numeric(), 
                     b0_coverage=numeric(), 
                     b0_type1=numeric(), 
                     b0_hat_bias=numeric(),
                     b1_coverage=numeric(), 
                     b1_type1=numeric(), 
                     b1_hat_bias=numeric())
for(beta0 in c(0)) {
  for(beta1 in c(0, 2.3)) {
    for(rho in c(0.9, 0.5, 0)) {
      print(paste0("beta0=",beta0," beta1=",beta1, " rho=", rho))
      for(m in c(10,50)) {
        for(n in c(50,10)) {
          if(m!=n) {  
            r <- run_sims(N=N, m=m, n=n, mu=0, rho=rho, beta1=beta1, beta0=0)
            tab <- dplyr::add_row(tab, 
                                  "beta0"=beta0, "beta1"=beta1, 
                                  "rho"=rho, "m"=m, "n"=n, 
                                  "b0_coverage"=r[[1]], 
                                  "b0_type1"=r[[3]], 
                                  "b0_hat_bias"=round(r[[5]],2),
                                  "b1_coverage"=r[[2]], 
                                  "b1_type1"=r[[4]], 
                                  "b1_hat_bias"=round(r[[6]],2)                                )
          }
        }
      }
    }
  }
}
tab

# ~~~~~~~~~ #
# Problem 2 #
# ~~~~~~~~~ #
#rm(list = ls())
#m=5; n=4; mu=0; rho=0.9; b1=1.9; alpha = 0.05

simul_logit <- function(m, n, mu=0, rho, b1, alpha=0.05) {
  
  # create correlated X which will determine probabilities p
  mus <- rep(mu, n)
  sigma <- matrix(rep(rho, n^2), nrow=n, ncol=n)
  diag(sigma) <- 1
  sigma
  b0 <- matrix(rep(runif(m,-1,1),n), nrow=m, ncol=n)
  X <- mvtnorm::rmvnorm(n=m, mean=mus, sigma=sigma)
  p <- exp(b0+b1*X)/(1+exp(b0+b1*X))
  Y <- rbinom(n=m*n, size=1, prob=p)
  id <- rep(1:m, n)
  df <- data.frame(id=id, Y=Y, 
                   X=matrix(X, nrow=m*n, ncol=1), 
                   p=matrix(p, nrow=m*n, ncol=1)) |> dplyr::arrange(id)
  df
  # fit logistic model
  cz <- qnorm(1-alpha/2, 0, 1)
  fit <- summary(glm(Y~X, data=df, family=binomial(link="logit")))$coefficients
  b1_hat    <- fit[2,1]
  inside    <- as.numeric(fit[2,1]-cz*fit[2,2] < b1 & b1 < fit[2,1]+cz*fit[2,2])
  signif_b1 <- as.numeric(fit[2,4] < 0.05)
  type_1    <- as.numeric(signif_b1==1 & b1==0)
  
  # output per iteration
  o <- list("inside" = inside, "type_1" = type_1, "b1_hat" = b1_hat)
  return(o)
}

run_logsims <- function(N, m, n, mu, rho, b1, alpha){
  coverage <- rep(0, N)
  type_1 <- rep(0,N)
  b1_hat <- rep(0,N)
  for(i in 1:N) {
    r <- simul_logit(m=m, n=n, mu=0, rho=rho, b1=b1, alpha=0.05)
    coverage[i] <- r[[1]] 
    type_1[i] <- r[[2]] 
    b1_hat[i] <- r[[3]]
  }
  o <- list("coverage" = mean(coverage), 
            "type_1" = mean(type_1), 
            "b1_hat" = round(mean(b1_hat),2))
  return(o)
}

tab_log <- dplyr::tibble(b1=numeric(), rho=numeric(), 
                         m=numeric(), n=numeric(), 
                         coverage=numeric(), 
                         TypeI=numeric(), 
                         Beta1_hat=numeric())
for(b1 in c(0, 1.9)) {
  for(rho in c(0.9, 0.5, 0)) {
    print(paste0("b1=",b1, " rho=", rho))
    for(m in c(10,50)) {
      for(n in c(50,10)) {
        if(m!=n) {  
          r <- run_logsims(N=N, m=m, n=n, mu=0, rho=rho, b1=b1, alpha=0.05)
          tab_log <- dplyr::add_row(tab_log, "b1"=b1, "rho"=rho, "m"=m, "n"=n, 
                                    "coverage"=r[[1]], 
                                    "TypeI"=r[[2]],
                                    "Beta1_hat"=r[[3]])
        }
      }
    }
  }
}
tab_log
tab_log$TypeI[tab_log$b1!=0] <-NA


objs <- ls()
objs <- ls()[!(ls() %in% c("tab","tab_log","simul_logit","simul_reg"))]
rm(list = objs)
tab
tab_log

#################### 
path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw1/hw1_objects.Rdata"))







# NOTES

# built in pipe ti create p in one fell swoop
#X <- mvtnorm::rmvnorm(n=1, mean=mus, sigma=sigma)
## use x to build probabilities
#p <- exp(b0+b1*X)/(1+exp(b0+b1*X))
# the same as: p <- mvtnorm::rmvnorm(n=1, mean=mus, sigma=sigma) |> (\(x) exp(b0+b1*x)/(1+exp(b0+b1*x)) )()
# see https://ivelasq.rbind.io/blog/understanding-the-r-pipe/

# Create correlated e vector with decaying correlations
### mu <- rep(mu, n)
### t <- rnorm(n, mean = 0, sd = 1)
### d <- expand.grid(t(t), t)
### covs <- rho^abs(d$Var1 - d$Var2)
### sigma <- matrix(covs, nrow=n, ncol=n)



# DEPRECATED CODE

# method 2
# # create correlated X which will determine probabilities p
# mus <- rep(mu, n)
# sigma <- matrix(rep(rho, n^2), nrow=n, ncol=n)
# diag(sigma) <- 1
# sigma
# X <- mvtnorm::rmvnorm(n=1, mean=mus, sigma=sigma)
# 
# # use x to build probabilities
# p <- exp(b0+b1*X)/(1+exp(b0+b1*X))
# 
# # use probabilities to spawn 0,1 values
# Y <- rbinom(n=n, size=1, prob=p)
# 
# # repeat for m individuals
# id <- rep(1, n)
# for(i in 2:m) {
#   mus <- rep(mu, n) # wiggle mean for each individual m
#   x <- mvtnorm::rmvnorm(n=1, mean=mus, sigma=sigma)
#   p <- exp(b0+b1*x)/(1+exp(b0+b1*x))
#   Y <- c(Y, rbinom(n=n, size=1, prob=p))
#   X <- c(X, x)
#   id <- c(id, rep(i, n))
# }
# data.frame(id=id, Y=Y, X=X) |> dplyr::mutate(p = exp(b0+b1*X)/(1+exp(b0+b1*X)))
#~~~~~~~~~~~~~~~~~~~~~

# method 1
# generate a MVN and then create 1,0 variable based on the MVN values and a probability of my choice
# how do I incorporate X and the true beta parameter???
# mu <- rep(mu, n)
# t <- rnorm(n, mean=0, sd=1)
# d <- expand.grid(t(t), t)
# covs <- rho^abs(d$Var1 - d$Var2)
# sigma <- matrix(covs, nrow=n, ncol=n)
# v1 <- mvtnorm::rmvnorm(n=1, mean=mu, sigma=sigma)
# y <- ifelse(v1>median(v1), 
#             sample(c(0,1),n, replace = TRUE, p=c(myp,1-myp)),
#             sample(c(0,1),n, replace = TRUE, p=c(1-myp,myp)))
# 
# 
# 
# # method 2
# 
# x <- mvtnorm::rmvnorm(n=m, mean=mu, sigma=sigma)
# p <- exp(b0+b1*x[1,])/(1+exp(b0+b1*x[1,]))
# Y <- rbinom(n=n, size=1, prob = exp(b0+b1*x[i,])/(1+exp(b0+b1*x[i,])))
# 
# Y <- y[1,]
# X <- x[1,]
# 
# for(i in 2:m) {
#   Y <- c(Y, y[i,]) #<- c(Y, rbinom(n=n, size=1, prob = exp(b0+b1*x[i,])/(1+exp(b0+b1*x[i,]))))
#   X <- c(X, x[i,])
# }
# 
# 
# cz <- qnorm(1-alpha/2, 0, 1)
# fit <- summary(glm(Y~X, family = binomial(link="logit")))$coefficients
# beta_1_hat <- fit[2,1]
# inside <- as.numeric(fit[2,1]-cz*fit[2,2] < b1 & b1 < fit[2,1]+cz*fit[2,2])
# signif_b1  <- as.numeric(fit[2,4] < 0.05)
# type_i <- as.numeric(signif_b1 == 1 & b1==0)
# 
# o <- list("inside" = inside, "signif_b1" = signif_b1)
