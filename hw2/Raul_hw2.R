library(mvtnorm)
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(lmerTest)
set.seed(5712)
rm(list = ls())
alpha=0.05
N <- 1e2


# ~~~~~~~~~ #
# Problem 2 #
# ~~~~~~~~~ #
m=200
n=3
beta0 = 1
beta1 = 0.5
theta = 1
sigma2 = 1

# create random intercept b_rand ~ N(0,theta)
b_r <- rep(rnorm(m, 0, theta),n)

# Create random noise
e <- rnorm(m*n, 0, sigma2)

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

# inspect dataset
df_p2 <- data.frame(id=id, Y=Y, X=X, b_r=b_r, e=e) |> dplyr::arrange(id)
for(i in 1:m) {
  df_p2[["temp"]] <- as.numeric(df_p2$id == i)
  colnames(df_p2)[ncol(df_p2)] <- paste0("Z",i)
}

mod_reml <- VarCorr(lme4::lmer(Y~X + (1|id), REML = TRUE)) |> as.data.frame()
mod_ml   <- VarCorr(lme4::lmer(Y~X + (1|id), REML = FALSE)) |> as.data.frame()

myReML_ll <- function(X, Y, Z, m, n, par) {
  theta <- par[1]
  sigma2 <- par[2]
  
  D <- theta * diag(m)
  R <- sigma2 * diag(m*n)
  V <- Z%*% D %*% t(Z) + R
  
  V_inv <- solve(V)
  beta <- solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv %*% Y
  W <- Y - X %*% beta
  
  l <- -0.5 * log(det(V)) -0.5*(t(W) %*% V_inv %*% W) -0.5 * log(det(t(X) %*% V_inv %*% X))
  return(-l)
}

myML_ll <- function(X, Y, Z, m, n, par) {
  theta <- par[1]
  sigma2 <- par[2]
  
  D <- theta * diag(m)
  R <- sigma2 * diag(m*n)
  V <- Z %*% D %*% t(Z) + R
  
  W <- Y - (X %*% solve((t(X)%*%X))%*%t(X)%*%Y)
  V_inv <- solve(V)
  l <- -0.5 * log(det(V)) -0.5 * (t(W) %*% V_inv %*% W)
  return(-l)
}


byhand_ReML <- optim(par = c(1,1),
                     fn = myReML_ll, 
                     X = as.matrix(cbind(1,df_p2$X)), 
                     Y = as.matrix(df_p2$Y), 
                     Z = as.matrix(df_p2 %>% dplyr::select(starts_with("Z"))), 
                     m=m, n=n)$par

byhand_ML   <- optim(par = c(1,1),
                     fn = myML_ll,
                     X = as.matrix(cbind(1,df_p2$X)), 
                     Y = as.matrix(df_p2$Y), 
                     Z = as.matrix(df_p2 %>% dplyr::select(starts_with("Z"))), 
                     m=m, n=n)$par

tabp2 <- dplyr::tibble(param = character(),
                       lmer_ML=numeric(), byhand_ML=numeric(), 
                       lmer_ReML=numeric(), byhand_ReML=numeric())
tabp2 <- dplyr::add_row(tabp2, "param"="sigma2", "lmer_ML"=mod_ml[1,4], "byhand_ML"=byhand_ML[2], 
                               "lmer_ReML"=mod_reml[1,4], "byhand_ReML"=byhand_ReML[2])
tabp2 <- dplyr::add_row(tabp2, "param"="theta" , "lmer_ML"=mod_ml[2,4], "byhand_ML"=byhand_ML[1], 
                               "lmer_ReML"=mod_reml[2,4], "byhand_ReML"=byhand_ReML[1])

plot(x=c(0.982, 0.979, 0.995, 0.979), y=c(1,2,3,4), 
     xlab = "estimate", 
     ylab = "",
     xlim = c(0.97, 1.02), 
     yaxt = "n",
     ylim = c(1,5))
abline(v = 1, col="red", lwd=3, lty=2)
     text(c(0.982, 0.979, 0.995, 0.979), 
          y=c(1,2,3,4)+.1, 
          labels=c("lmer ML","by hand ML","lmer ReML","by hand ReML"))


# ~~~~~~~~~ #
# Problem 3 #
# ~~~~~~~~~ #

simul_reg <- function(m, n, rho, beta0=1, beta1=0, sigma=1) {
  
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
  
  # Fit linear mixed model (ML)
  fit <- lme4::lmer(Y~X + (1|id), REML = FALSE)
  beta_hat <- coef(summary(fit))
  CI_beta_hat_ml <- c(beta_hat[2,1] - beta_hat[2,3], beta_hat[2,1] + beta_hat[2,3])
  inside_ml <- as.numeric(CI_beta_hat_ml[1] <= beta1 & beta1 <= CI_beta_hat_ml[2])
  var_hat_ml <- as.data.frame(VarCorr(fit))[1,5]
  
  # Fit linear mixed model (ML)
  fit_reml <- lme4::lmer(Y~X + (1|id), REML = TRUE)
  beta_hat_re <- coef(summary(fit_reml))
  CI_beta_hat_reml <- c(beta_hat_re[2,1] - beta_hat_re[2,3], beta_hat_re[2,1] + beta_hat_re[2,3])
  inside_reml <- as.numeric(CI_beta_hat_reml[1] <= beta1 & beta1 <= CI_beta_hat_reml[2])
  var_hat_reml <- as.data.frame(VarCorr(fit_reml))[1,5]
  
  o <- list("true_var" = theta,
            "b1_inside_ML" = inside_ml[[1]], 
            "CI_ML_025" = CI_beta_hat_ml[[1]], 
            "CI_ML_975" = CI_beta_hat_ml[[2]],
            "var_ML"    = var_hat_ml[[1]], 
            "b1_inside_ReML" = inside_reml[[1]], 
            "CI_ReML_025" = CI_beta_hat_reml[[1]], 
            "CI_ReML_975" = CI_beta_hat_reml[[2]],
            "var_ReML"    = var_hat_reml[[1]])
  return(o)
}

simul_reg(m=5, n=3, rho=0.85, beta0=1, beta1=2, sigma=1)

# ~~~~~~~~~~~~~~~~~~~~~~~ #
# create table of results #
# ~~~~~~~~~~~~~~~~~~~~~~~ #

run_sims <- function(N, m, n, rho, beta0, beta1, alpha=0.05, sigma=1) {
  true_var     <- rep(0,N)
  b1_inside_ML <- rep(0,N) ;b1_inside_ReML <- rep(0,N)
  CI_ML_025    <- rep(0,N) ;CI_ReML_025    <- rep(0,N); var_ML   <- rep(0,N)
  CI_ML_975    <- rep(0,N) ;CI_ReML_975    <- rep(0,N); var_ReML <- rep(0,N)
  for(i in 1:N) {
    r <- simul_reg(m=m, n=n, rho=rho, beta0=beta0, beta1=beta1, sigma=sigma)
    true_var[i]         <- r[[1]] 
    b1_inside_ML[i]     <- r[[2]]
    CI_ML_025[i]        <- r[[3]] 
    CI_ML_975[i]        <- r[[4]]
    var_ML[i]           <- r[[5]] 
    b1_inside_ReML[i]   <- r[[6]] 
    CI_ReML_025[i]      <- r[[7]]    
    CI_ReML_975[i]      <- r[[8]]
    var_ReML[i]         <- r[[9]]
  }
  o <- list("true_var"       = mean(true_var),
            "b1_inside_ML"   = mean(b1_inside_ML),
            "CI_ML_025"      = mean(CI_ML_025),
            "CI_ML_975"      = mean(CI_ML_975),
            "var_ML"         = mean(var_ML)-true_var[1],
            "b1_inside_ReML" = mean(b1_inside_ReML),
            "CI_ReML_025"    = mean(CI_ReML_025),
            "CI_ReML_975"    = mean(CI_ReML_975),
            "var_ReML"       = mean(var_ReML)-true_var[1]            
            )
  return(o)
}

run_sims(1e2, m=5,n=3, rho=0.85, beta0=1, beta1=2, sigma=1)

tab <- dplyr::tibble(rho=numeric(), m=numeric(), n=numeric(), N=numeric(),
                     variance=numeric(),
                     bias_var_ML=numeric(), 
                     bias_var_ReML=numeric())
for(rho in c(0.85)) {
  for(m in c(5, 20, 50)) {
    for(n in c(3, 10, 30)) {
      r <- run_sims(N=1e2, m=m, n=n, rho=rho, beta1=1.7, beta0=1, sigma=1)
      tab <- dplyr::add_row(tab,
                            "rho"=rho, "m"=m, "n"=n, "N"=m*n,
                            "variance" = r[[1]],
                            "bias_var_ML" = r[[5]],
                            "bias_var_ReML" = r[[9]])
    }
  }
}



# ~~~~~~~~~ #
# Problem 4 #
# ~~~~~~~~~ #

df <- read.table("hw2/framingham.dat") |> 
  `colnames<-`(c("age", "gender", "bmi_0", "bmi_10", "cigs", 
                 "col_0","col_2","col_4","col_6","col_8","col_10","dead")) |>
  mutate_if(is.numeric, ~replace(., . == -9, NA)) |>
  mutate(id = row_number()) |> 
  dplyr::select(id, dead, everything()) |>
  pivot_longer(cols = starts_with("col_"),
               names_to = "t",
               names_pattern = "col_(.+)",
               values_to = "cholesterol") 

df$bmi <- NA
df$bmi[df$t==0] <- df$bmi_0[df$t==0] 
df$bmi[df$t==10] <- df$bmi_10[df$t==10] 
df$t <- as.numeric(df$t)
df <- df |> 
  dplyr::select(id, t, dead, cholesterol, bmi, bmi_0, bmi_10, everything())
df$gender <- factor(df$gender)


  # naive plot all
  plot(df$cholesterol~df$t, xlab="time", ylab="cholesterol", xaxt="n")
  abline(lm(df$cholesterol~df$t), col='red')
  axis(1, at = seq(from=0, to=10, by=2), 
       labels = c("baseline","2y","4y","6y","8y","10y")) 
  
  # naive plot all
  plot(df$cholesterol~df$t, xlab="time", ylab="cholesterol", xaxt="n")
  abline(lm(df[df$gender == 1,]$cholesterol~df[df$gender == 1,]$t), col='lightblue', lwd=4.0)
  abline(lm(df[df$gender == 2,]$cholesterol~df[df$gender == 2,]$t), col='pink',lwd=4.0)
  axis(1, at = seq(from=0, to=10, by=2), 
       labels = c("baseline","2y","4y","6y","8y","10y")) 
  

  # naive plot 5 randomly selected
  rand_id <- sample(1:max(df$id), 50, replace = FALSE)
  df_rand_all <- df[df$id %in% rand_id, ]
  plot(df_rand$cholesterol~df_rand$t, xlab="time", ylab="cholesterol", xaxt="n", main = "five random cases")
  axis(1, at = seq(from=0, to=10, by=2), 
       labels = c("baseline","2y","4y","6y","8y","10y"))
  abline(lm(df_rand_all$cholesterol[df_rand_all$id == rand_id[1]]~df_rand_all$t[df_rand$id == rand_id[1]]), col='blue')
  abline(lm(df_rand_all$cholesterol[df_rand_all$id == rand_id[2]]~df_rand_all$t[df_rand$id == rand_id[2]]), col='orange')
  abline(lm(df_rand_all$cholesterol[df_rand_all$id == rand_id[3]]~df_rand_all$t[df_rand$id == rand_id[3]]), col='red')
  abline(lm(df_rand_all$cholesterol[df_rand_all$id == rand_id[4]]~df_rand_all$t[df_rand$id == rand_id[4]]), col='darkgreen')
  abline(lm(df_rand_all$cholesterol[df_rand_all$id == rand_id[6]]~df_rand_all$t[df_rand$id == rand_id[6]]), col='purple')

  
  par(mfrow=c(1,2))
  # age and sex and cholesterol
  plot(df$cholesterol~df$age, xlab="age", ylab="cholesterol", xaxt="n")
  abline(lm(df[df$gender == 1,]$cholesterol~df[df$gender == 1,]$age), col='lightblue', lwd=4.0)
  abline(lm(df[df$gender == 2,]$cholesterol~df[df$gender == 2,]$age), col='pink',lwd=4.0)

  # age and sex and cholesterol random sample
  df_rand <- df_rand_all
  df_rand$age <- df_rand$age + df_rand$t
  df_rand <- df_rand[!is.na(df_rand$cholesterol),]
  df_rand <- df_rand %>% group_by(id) %>% mutate(size = n()) %>% filter(size>5)
  example_male   <- unique(df_rand$id[df_rand$gender==1])[1:10]
  example_female <- unique(df_rand$id[df_rand$gender==2])[1:10]
  plot(df_rand$cholesterol~df_rand$age, xlab="age", ylab="cholesterol", main = "chol over age by sex:\n6 random subjects")
  abline(lm(df_rand[df_rand$gender == 1 & df_rand$id == example_male[1],]$cholesterol~
            df_rand[df_rand$gender == 1 & df_rand$id == example_male[1],]$age), col='lightblue', lwd=4.0)
  abline(lm(df_rand[df_rand$gender == 2 & df_rand$id == example_female[1],]$cholesterol~
            df_rand[df_rand$gender == 2 & df_rand$id == example_female[1],]$age), col='pink',lwd=4.0)
  
  abline(lm(df_rand[df_rand$gender == 1 & df_rand$id   == example_male[2]  ,]$cholesterol~
              df_rand[df_rand$gender == 1 & df_rand$id == example_male[2]  ,]$age), col='lightblue', lwd=4.0)
  abline(lm(df_rand[df_rand$gender == 2 & df_rand$id   == example_female[2],]$cholesterol~
              df_rand[df_rand$gender == 2 & df_rand$id == example_female[2],]$age), col='pink',lwd=4.0)
  
  abline(lm(df_rand[df_rand$gender == 1   & df_rand$id == example_male[3],]$cholesterol~
              df_rand[df_rand$gender == 1 & df_rand$id == example_male[3],]$age), col='lightblue', lwd=4.0)
  abline(lm(df_rand[df_rand$gender == 2.  & df_rand$id == example_female[3],]$cholesterol~
              df_rand[df_rand$gender == 2 & df_rand$id == example_female[3],]$age), col='pink',lwd=4.0)    
  par(mfrow=c(1,1))
  

# (4.1)
# how does cholesterol change over time and how is it related to age, gender, BMI?
#fit1 <- lme4::lmer(cholesterol ~ t + t*age + t*gender + t*bmi_0 + t*bmi_10 + (1+t|id), REML=FALSE, data=df)
fit1 <- lme4::lmer(cholesterol ~ t + age*gender + bmi_0 + bmi_10 + (1+t|id), REML=FALSE, data=df)

sfit1   <- summary(fit1)
coeffs1 <- sfit1$coefficients #|> round(2)
corrs1  <- cov2cor(vcov(fit1)) #|> round(2)
sfit1

# retrieve random effect of t per id
REs <- data.frame(id = 1:length(coef(fit1)$id$t), 
                  re_t = coef(fit1)$id$t, 
                  re_id = coef(fit1)$id$`(Intercept)`)
df2 <- df |> 
  left_join(REs, by = "id") |> 
  group_by(id) |>
  mutate(chol_0 = first(cholesterol))
print(df2, n=18)

# (4.2)
fit2 <- glm(dead ~ re_id + re_t + age*gender + bmi_0 + bmi_10 + chol_0, 
            data = df2[df2$t==0,],
            family = "binomial") 
sfit2   <- summary(fit2)
coeffs2 <- sfit2$coefficients 
coeffs2[,1] <- exp(coeffs2[,1]) 
coeffs2 |> round(3)

# (4.3)


#################### 
path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw2/hw2_objects.Rdata"))







