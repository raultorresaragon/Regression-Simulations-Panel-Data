library(mvtnorm)
set.seed(571)
rm(list = ls())
alpha=0.05
N <- 1e3

# ~~~~~~~~~ #
# Problem 1 #
# ~~~~~~~~~ #
#m=5; n=4; mu=0; rho=.78; b0=1; b1=2.5

simul_reg <- function(m, n, mu=0, rho, b0=1, b1=0, alpha=0.05) {
  
  # Create correlated e vector with constant correlations
  mu <- rep(mu, n)
  Sigma <- matrix(rep(rho, n^2), nrow=n, ncol=n)
  diag(sigma) <- 1
  Sigma
  e <- mvtnorm::rmvnorm(n=m, mean=rep(mu,n), sigma=Sigma) |>
       matrix(ncol=1, nrow=m*n)
  id <- rep(1:m,n)
  
  
  e <- mvtnorm::rmvnorm(n=1, mean=mu, sigma=sigma)
  id <- rep(1, n)
  for(i in 2:m) {
    e <- c(e, mvtnorm::rmvnorm(n=1, mean=mu, sigma=sigma))
    id <- c(id, rep(i, n))
  }
  
  # Create X vector
  trt <- sample(c(0,1),1)
  X <- rep(trt, m*n)
  if(trt==1) { 
    X[id %% 2 == 0] <- 0 
  }else { 
    X[id %% 2 == 0] <- 1 
  }

  # Create Y vector with true relationship
  Y <- b0 + b1*X + e 
  
  # Fit linear model
  fit <- lm(Y~X)
  Beta <- coef(fit)
  Beta_SE <- sqrt(diag(vcov(fit)))
  CI <- c(Beta[2] - abs(qt(alpha/2, m*n-2))*Beta_SE[2], 
          Beta[2] + abs(qt(alpha/2, m*n-2))*Beta_SE[2])
  #CI <- confint(fit)[2,]
  inside <- as.numeric(CI[1] <= b1 & b1 <= CI[2])
  if(verbose==TRUE) { print(CI) }
  
  pval <- summary(fit)$coefficients[2,4]
  signific <- as.numeric(pval<0.05)
  if(verbose==TRUE) { print(pval) }

  #type I error: rejecting the null when the null is true
  if(b1==0) {
    type_i <- as.numeric(signific==1 & b1==0)
  } else{
    type_i <- as.numeric(signific==0 & b1!=0)
  }
    
  o <- list("inside" = inside, "type_i" = type_i, "b1hat" = Beta[2])
  return(o)
}


#################### create table of results
run_sims <- function(N, m, n, mu, rho, b0, b1, alpha=0.05) {
  coverage <- rep(0, N)
  type_i   <- rep(0,N)
  b1_hat    <- rep(0,N)
  for(i in 1:N) {
    r <- simul_reg(m, n, mu, rho, b0, b1, alpha)
    coverage[i] <- r[[1]] 
    type_i[i] <- r[[2]] 
    b1_hat[i] <- r[[3]]
  }
  o <- list("coverage" = mean(coverage), "type_i" = mean(type_i), "b1_hat" = mean(b1_hat))
  return(o)
}

tab <- dplyr::tibble(b1=numeric(), rho=numeric(), 
                     m=numeric(), n=numeric(), 
                     coverage=numeric(), TypeI=numeric(), Beta1_hat=numeric())
for(b1 in c(0, 2.3)) {
  for(rho in c(0.9, 0.5, 0)) {
    print(paste0("b1=",b1, " rho=", rho))
    for(m in c(10,50)) {
      for(n in c(50,10)) {
        if(m!=n) {  
        r <- run_sims(N=N, m=m, n=n, mu=0, rho=rho, b1=b1, b0=0, alpha=0.05)
        tab <- dplyr::add_row(tab, "b1"=b1, "rho"=rho, "m"=m, "n"=n, 
                              "coverage"=r[[1]], 
                              "TypeI"=r[[2]], 
                              "Beta1_hat"=round(r[[3]],2))
        }
      }
    }
  }
}
tab$TypeI[tab$b1!=0] <-NA
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
