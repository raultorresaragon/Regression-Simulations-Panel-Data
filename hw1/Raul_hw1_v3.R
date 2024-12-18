library(mvtnorm)
set.seed(571)
rm(list = ls())
alpha=0.05
N <- 1e3

# ~~~~~~~~~ #
# Problem 1 #
# ~~~~~~~~~ #

simul_reg <- function(m, n, rho, beta0=1, beta1=0, sigma=1, alpha=0.05) {
  
  # create random intercept b_rand ~ N(0,theta)
  theta = -(rho*sigma^2)/(rho-1)
  b_r <- rep(rnorm(m, 0, theta),n)
  
  # Create correlated e vector
  e <- rnorm(m*n, 0, sigma)
  
  # Create id vector for sanity checks
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
  Y <- beta0 + beta1*X + b_r + e 
  
  data.frame(id=id, Y=Y, X=X, b_r=b_r, e=e) |> dplyr::arrange(id)
  
  # Fit linear model
  fit <- lm(Y~X)
  beta_hat <- coef(fit)
  beta_hat_SE <- sqrt(diag(vcov(fit)))
  CI_beta_hat <- matrix(c(beta_hat - abs(qt(alpha/2, m*n-2))*beta_hat_SE, 
                          beta_hat + abs(qt(alpha/2, m*n-2))*beta_hat_SE), ncol=2, nrow=2)
  #CI_beta_hat <- confint(fit)
  inside <- c(as.numeric(CI_beta_hat[1,1] <= beta0 & beta0 <= CI_beta_hat[1,2]),
              as.numeric(CI_beta_hat[2,1] <= beta1 & beta1 <= CI_beta_hat[2,2]))
  pval <- summary(fit)$coefficients[,4]
  
  #type 1 error: rejecting the null when the null is true
  type_1 <- c(NA,NA)
  if(beta0==0) {
    type_1[1] <- as.numeric(pval[1]<0.05)
  }
  if(beta1==0) {
    type_1[2] <- as.numeric(pval[2]<0.05)
  }
  
  o <- list("b0_inside" = inside[[1]], 
            "b0_type_1" = type_1[[1]], 
            "b0_hat"    = beta_hat[[1]],
            "b1_inside" = inside[[2]], 
            "b1_type_1" = type_1[[2]], 
            "b1_hat"    = beta_hat[[2]])
  return(o)
}

#~~~~~~~~~~~#
# Problem 2 #
#~~~~~~~~~~~#
simul_logit <- function(m, n, rho, beta0=1, beta1=0, alpha=0.05, sigma=1) {
  
  # create correlated X which will determine probabilities p
  theta = -(rho*sigma^2)/(rho-1)
  b_r <- rep(rnorm(m, 0, theta),n)
  X <- rnorm(n*m, 0,1)
  p <- exp(beta0+b_r+beta1*X)/(1+exp(beta0+b_r+beta1*X))
  Y <- rbinom(n=m*n, size=1, prob=p)
  id <- rep(1:m, n)
  df <- data.frame(id=id, Y=Y, 
                   X=matrix(X, nrow=m*n, ncol=1), 
                   p=matrix(p, nrow=m*n, ncol=1)) |> dplyr::arrange(id)
  df

  # fit logistic model
  cz <- qnorm(1-alpha/2, 0, 1)
  fit <- summary(glm(Y~X, data=df, family=binomial(link="logit")))$coefficients
  beta_hat  <- fit[,1]
  inside    <- c(as.numeric(fit[1,1]-cz*fit[1,2] < beta0 & beta0 < fit[1,1]+cz*fit[1,2]),
                 as.numeric(fit[2,1]-cz*fit[2,2] < beta1 & beta1 < fit[2,1]+cz*fit[2,2]))
  pval <- fit[,4]
  
  #type 1 error: rejecting the null when the null is true
  type_1 <- c(NA,NA)
  if(beta0==0) {
    type_1[1] <- as.numeric(pval[1]<0.05)
  }
  if(beta1==0) {
    type_1[2] <- as.numeric(pval[2]<0.05)
  }
  
  # output per iteration
  o <- list("b0_inside" = inside[[1]], 
            "b0_type_1" = type_1[[1]], 
            "b0_hat"    = beta_hat[[1]],
            "b1_inside" = inside[[2]], 
            "b1_type_1" = type_1[[2]], 
            "b1_hat"    = beta_hat[[2]])
  return(o)
}



# ~~~~~~~~~~~~~~~~~~~~~~~ #
# create table of results #
# ~~~~~~~~~~~~~~~~~~~~~~~ #

run_sims <- function(N, m, n, rho, beta0, beta1, alpha=0.05, sigma=1, flavor="ols") {
  coverage_b0 <- rep(0,N)    ;coverage_b1 <- rep(0,N)
  type_1_b0   <- rep(0,N)    ;type_1_b1   <- rep(0,N)
  b0_hat      <- rep(0,N)    ;b1_hat      <- rep(0,N)
  for(i in 1:N) {
    if(flavor=="logit") {
      r <- simul_logit(m=m, n=n, rho=rho, beta1=beta1, beta0=beta0, sigma=sigma)
    } else{
      r <-   simul_reg(m=m, n=n, rho=rho, beta0=beta0, beta1=beta1, sigma=sigma)
    }
    coverage_b0[i] <- r[[1]] 
    type_1_b0[i]   <- r[[2]] 
    b0_hat[i]      <- r[[3]]
    coverage_b1[i] <- r[[4]] 
    type_1_b1[i]   <- r[[5]] 
    b1_hat[i]      <- r[[6]]    
  }
  o <- list("coverage_b0" = mean(coverage_b0) ,"coverage_b1" = mean(coverage_b1),
            "type_1_b0"   = mean(type_1_b0)   ,"type_1_b1"   = mean(type_1_b1),
            "b0_hat"      = beta0-mean(b0_hat),"b1_hat"      = beta1-mean(b1_hat))
  return(o)
}

tab <- dplyr::tibble(beta0=numeric(), beta1=numeric(),
                     rho=numeric(), m=numeric(), n=numeric(), 
                     b0_coverage=numeric(), 
                     b0_type1=numeric(), 
                     b0_hat_bias=numeric(),
                     b1_coverage=numeric(), 
                     b1_type1=numeric(), 
                     b1_hat_bias=numeric())

tab_logit <- dplyr::tibble(beta0=numeric(), beta1=numeric(),
                     rho=numeric(), m=numeric(), n=numeric(), 
                     b0_coverage=numeric(), 
                     b0_type1=numeric(), 
                     b0_hat_bias=numeric(),
                     b1_coverage=numeric(), 
                     b1_type1=numeric(), 
                     b1_hat_bias=numeric())

for(beta0 in c(0,1.9)) {
  for(beta1 in c(2.3)) {
    for(rho in c(0.9, 0.5, 0)) {
      print(paste0("beta0=",beta0," beta1=",beta1, " rho=", rho))
      for(m in c(10,50)) {
        for(n in c(50,10)) {
          if(m!=n) {  
            r <- run_sims(N=N, m=m, n=n, rho=rho, beta1=beta1, beta0=beta0, flavor="ols")
            tab <- dplyr::add_row(tab, 
                                  "beta0"=beta0, "beta1"=beta1, 
                                  "rho"=rho, "m"=m, "n"=n, 
                                  "b0_coverage"=r[[1]], 
                                  "b0_type1"=r[[3]], 
                                  "b0_hat_bias"=round(r[[5]],2),
                                  "b1_coverage"=r[[2]], 
                                  "b1_type1"=r[[4]], 
                                  "b1_hat_bias"=round(r[[6]],2))
            r <- run_sims(N=N, m=m, n=n, rho=rho, beta1=beta1, beta0=beta0, flavor="logit")
            tab_logit <- dplyr::add_row(tab_logit, 
                                  "beta0"=beta0, "beta1"=beta1, 
                                  "rho"=rho, "m"=m, "n"=n, 
                                  "b0_coverage"=r[[1]], 
                                  "b0_type1"=r[[3]], 
                                  "b0_hat_bias"=round(r[[5]],2),
                                  "b1_coverage"=r[[2]], 
                                  "b1_type1"=r[[4]], 
                                  "b1_hat_bias"=round(r[[6]],2))            
            
          }
        }
      }
    }
  }
}
tab
tab_logit
objs <- ls()
objs <- ls()[!(ls() %in% c("tab","tab_logit","simul_logit","simul_reg","run_sims"))]
rm(list = objs)


#################### 
path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw1/hw1_objects.Rdata"))







