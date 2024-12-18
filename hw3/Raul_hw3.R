library(dplyr)
library(tidyr)
library(gee)
library(lme4)
library(SimCorMultRes)
set.seed(571)
rm(list = ls())
df <- readr::read_table("hw3/sixcity.dat", col_names = c("ws","id","age","smokes"))

# ~~~~~~~~~ #
# Problem 1 #
# ~~~~~~~~~ #
myformula <- as.formula("ws ~ age +smokes +age:smokes")

### (1)  
# Fit a logistic reg ignoring within-subject correlation and calculate the naive SEs.  
mod1 <- glm(myformula, data=df, family=binomial(link="logit"))
coef1 <- as.data.frame(summary(mod1)$coefficients) 
coef1 

### (2)  
# Analyze this data set using GEE1 assuming working independence and exchangeable. 
# Compare the results with the naive logistic regression. 
# Interpret each regression coefficient, except for the intercept.
mod2 <- gee(myformula, data=df, id=id, family=binomial, corstr="independence")
coef2 <- as.data.frame(summary(mod2)$coefficients) 
coef2

mod3 <- gee(myformula, data=df, id=id, family=binomial, corstr="exchangeable")
coef3 <- as.data.frame(summary(mod3)$coefficients) 
coef3

tab1 <- cbind(coef1[,c(1,2)], coef2[,c(2,4)], coef3[c(4)])
colnames(tab1) <- c("Estimate","GLM","Naive","Robust (ind)","Robust (ex)")
tab1

# working plot creating confidence intervals: (but should we use n or N?)
# tab1$Estimate[1] - qnorm(0.975)*(tab1$`GLM SE`[1]/sqrt(nrow(df)))
# tab1$Estimate[1] + qnorm(0.975)*(tab1$`GLM SE`[1]/sqrt(nrow(df)))
#plot(tab1$Estimate[1], 1, pch=19, xlim = c(tab1$Estimate[1]-2,tab1$Estimate[1]+2))
#points(tab1$`GLM SE`[1], 1.1, pch=4, col="gray")
#points(tab1$`Naive SE`[1], 1.2, pch=4, col="darkgreen")
#points(tab1$`Robust_ind SE`[1], 1.3, pch=4, col="darkorange")
#points(tab1$`Robust_ex SE`[1], 1.4, pch=4, col="darkred")




### (3) (BONUS)  
#Analyze this data set using alternating logistic regression. Is there a strong within-subject correlation?

  
#### (4)  
#Propose a GEE1 model that can separate the baseline age effect from the longitudinal time effect.
df <- df |> group_by(id) |> mutate(time = row_number())
df$age <- df$age+2

#
mod4 <- gee(ws ~ factor(age) + smokes +smokes:age, corstr="exchangeable", id=id, family=binomial, data=df) 

# \text{ws}_{it} = \beta_0 + \beta_1 \times \text{age}_{base} + 
#                            \beta_2 \text{age} + \beta_3 text{smokes} + 
#                            \beta_4 \text{age} \times \text{smokes}
mod4 <- gee(ws ~ (age==0) + smokes + smokes*age, id=id, corstr="exchangeable", family=binomial, data=df,)

coef4 <- as.data.frame(summary(mod4)$coefficients) 
coef4
colnames(coef4) <- c("Estimate","Naive SE", "Naive Z", "Robust (ex) SE", "Robust Z")
coef4$OddsRatio <- exp(coef4$Estimate)



# ~~~~~~~~~ #
# Problem 2 #
# ~~~~~~~~~ #
m = 300
n = 3
rho = 0.25
beta0 = -1.5
beta1 = 0.5
beta2 = 0.5

## (2.1)
sim <- function(m=300, n=3, rho, beta0=-1.5, beta1=0.5, beta2=0.5) {
  
  cormatrix <- toeplitz(c(1,rho,rho))
  
  # id variable
  id <- rep(1:m, each=n)
  
  # Create X1 vector
  sex <- sample(c(0,1),1)
  X1 <- rep(sex, m*n)
  if(sex==1) { 
     X1[id %% 2 == 0] <- 0 
   }else { 
     X1[id %% 2 == 0] <- 1 
  }
  
  # Create X2 vector
  X2 <- rep(1:n, m)
  
  # Create Y vector
  Y <- rbin(clsize = n, intercepts = beta0, betas = c(beta1, beta2), 
            xformula = ~X1+X2, cor.matrix = cormatrix, link = "logit")$simdata$y
  
  # Pack them all in one dataframe
  df <- data.frame(id=id, Y=Y, X1=X1, X2=X2) |> arrange(id)
  
  # Run GEE
  mod    <- gee(Y ~ X1 + X2, corstr="exchangeable", id=id, data=df, family="binomial") 
  modsum <- summary(mod)
  modrho <- modsum$working.correlation[1,2]
  mod_b0 <- modsum$coefficients[1,1]
  mod_se0<- modsum$coefficients[1,4]
  mod_b1 <- modsum$coefficients[2,1]
  mod_se1<- modsum$coefficients[2,4]
  mod_b2 <- modsum$coefficients[3,1]
  mod_se2<- modsum$coefficients[3,4]
  
  o <- list("rho" = modrho, 
            "b0"  = mod_b0,
            "se0" = mod_se0,
            "b1"  = mod_b1,
            "se1" = mod_se1,
            "b2"  = mod_b2,
            "se2" = mod_se2)
  return(o)
}


  ## run 200 times
tab <- tibble(rho = numeric(),
              b0  = numeric(),
              b1  = numeric(),
              b2  = numeric(),
              se0 = numeric(),
              se1 = numeric(),
              se2 = numeric())
for(i in 1:200) {
  mysim <- sim(m=300, n=3, rho=0.25, beta0=-1.5, beta1=0.5, beta2=0.5)
  tab <- add_row(tab, 
                 "rho" = mysim$rho, 
                 "b0"  = mysim$b0,
                 "b1"  = mysim$b1,
                 "b2"  = mysim$b2,
                 "se0" = mysim$se0,
                 "se1" = mysim$se1,
                 "se2" = mysim$se2)
}
apply(tab, 2, mean)

## (2.2)
bias_tab <- data.frame(param = c(rho, beta0, beta1, beta2),
                       estimate = c(mean(tab$rho), mean(tab$b0), mean(tab$b1), mean(tab$b2)))
bias_tab$bias <- bias_tab$param - bias_tab$estimate
row.names(bias_tab) <- c("rho","beta0","beta1","beta2")
bias_tab


## (2.3) 
SE_tab <- data.frame(
            "Sandwich" = c(mean(tab$se0),mean(tab$se1),mean(tab$se2)),
            "Empirical"= c(sd(tab$b0), sd(tab$b1), sd(tab$b2))
          )
row.names(SE_tab) <- c("std. error beta0","std. error beta1","std. error beta2")

## (2.4)  
tab_24 <- tibble(rho = numeric(),
              b0  = numeric(),
              b1  = numeric(),
              b2  = numeric(),
              se0 = numeric(),
              se1 = numeric(),
              se2 = numeric())
for(i in 1:200) {
  mysim <- sim(m=300, n=3, rho=0.75, beta0=-1.5, beta1=0.5, beta2=0.5)
  tab_24 <- add_row(tab, 
                 "rho" = mysim$rho, 
                 "b0"  = mysim$b0,
                 "b1"  = mysim$b1,
                 "b2"  = mysim$b2,
                 "se0" = mysim$se0,
                 "se1" = mysim$se1,
                 "se2" = mysim$se2)
}
apply(tab_24, 2, mean)

bias_tab_24 <- data.frame(param = c(0.75, beta0, beta1, beta2),
                       estimate = c(mean(tab_24$rho), mean(tab_24$b0), mean(tab_24$b1), mean(tab_24$b2)))
bias_tab_24$bias <- bias_tab_24$param - bias_tab_24$estimate
row.names(bias_tab_24) <- c("rho","beta0","beta1","beta2")
bias_tab_24





# ~~~~~~~~~ #
# Problem 3 #
# ~~~~~~~~~ #

## (3.1)  
df_fram <- read.table("hw2/framingham.dat") |> 
           `colnames<-`(c("age", "gender", "bmi_0", "bmi_10", "cigs", 
                          "col_0","col_2","col_4","col_6","col_8","col_10","dead")) |>
           mutate_if(is.numeric, ~replace(., . == -9, NA)) |>
           mutate(id = row_number()) |> 
           dplyr::select(id, dead, everything()) |>
           pivot_longer(cols = starts_with("col_"),
                        names_to = "t",
                        names_pattern = "col_(.+)",
                        values_to = "cholesterol") |> 
           mutate(t = as.numeric(t),
                  gender = factor(gender),
                  age = age + t)

df_fram$bmi <- NA
df_fram$bmi[df_fram$t==0] <- df_fram$bmi_0[df_fram$t==0] 
df_fram$bmi[df_fram$t==10] <- df_fram$bmi_10[df_fram$t==10] 
df_fram <- df_fram |> 
           dplyr::select(id, t, dead, cholesterol, bmi, bmi_0, bmi_10, everything())


  ## LMM
  # (from hw2)
  # how does cholesterol change over time and how is it related to age, gender, BMI?
  fit_hw2 <- lme4::lmer(cholesterol ~ t + age*gender + bmi_0 + bmi_10 + (1+t|id), REML=FALSE, data=df_fram)
  sfit_hw2   <- summary(fit_hw2)
  coeffs1 <- sfit_hw2$coefficients #|> round(2)
  corrs1  <- cov2cor(vcov(fit_hw2)) #|> round(2)
  sfit_hw2
  # retrieve random effect of t per id
  REs <- data.frame(id = 1:length(coef(fit_hw2)$id$t), 
                    re_t = coef(fit_hw2)$id$t, 
                    re_id = coef(fit_hw2)$id$`(Intercept)`)
  df2 <- df_fram |> 
         left_join(REs, by = "id") |> 
         group_by(id) |>
         mutate(chol_0 = first(cholesterol))
  print(df2, n=18)

  # fit glm
  fit2 <- glm(dead ~ re_id + re_t + age*gender + bmi_0 + bmi_10 + chol_0, 
              data = df2[df2$t==0,],
              family = "binomial") 
  sfit2   <- summary(fit2)
  coeffs_hw2 <- sfit2$coefficients 
  coeffs_hw2[,1] <- exp(coeffs_hw2[,1]) 
  coeffs_hw2 |> round(3)

  ## GEE
  fit_gee <- gee(dead ~ t + age*gender + bmi_0 + bmi_10 + chol_0, 
                 corstr="exchangeable", 
                 id=id, 
                 family="binomial",
                 data=df2)
  fit_gee_coeffs <- summary(fit_gee)$coefficients[,c(1,4,5)]
  
  ## GEE ind
  fit_gee_ind <- gee(dead ~ t + age*gender + bmi_0 + bmi_10 + chol_0, 
                 corstr="independence", 
                 id=id, 
                 family="binomial",
                 data=df2)
  fit_gee_coeffs_ind <- summary(fit_gee_ind)$coefficients[,c(1,4,5)]
  
  
  par(mar=c(8,8,1,1)) 
  plot(exp(fit_gee_coeffs[2:8,1]), y=1:7, pch=1, col="darkblue", yaxt="n", ylab="", xlab="odds ratio") 
  points(x=coeffs_hw2[3:9,1]       , y=1:7, pch=4, col="darkred")
  points(exp(fit_gee_coeffs_ind[2:8,1]), y=1:7, pch=5, col="darkgreen") 
  axis(side=2, 
       at=1:7, 
       labels = c("t (re_t)", "age", "female", "bmi_baseline", "bmi_10wks", "chol_baseline", "age:female"),
       las=1
  )
  abline(v=1, col="black")
  grid(nx = NA, ny = NULL,
       lty = 2,      # Grid line type
       col = "gray", # Grid line color
       lwd = 2)  
  legend("topleft", legend=c("LME+logistic", "GEE (ex)", "GEE (in)"),
         col=c("darkblue", "darkred","darkgreen"), pch=c(1,4,5), cex=0.8)
  
  tab_31 <- data.frame("GEE_in" = exp(fit_gee_coeffs_ind[2:8,1]),
                       "GEE_ex" = exp(fit_gee_coeffs[2:8,1]),
                       "LMM" = coeffs_hw2[3:9,1]) |> mutate(diff = GEE_in - LMM)
                       
  row.names(tab_31) <- c("t (re_t)", "age", "female", "bmi (baseline)", "bmi (10 weeks)", "chol (baseline)", "age x female")
  colnames(tab_31) <- c("GEE (in)", "GEE (ex)", "LMM", "Diff (ex)")

## (3.2)  
sim_p32 <- function(m, n, beta0, beta1, rho) {
    
    # ID vector
    id <- rep(1:m, each=n)
    
    # Create X1 vector
    trt <- sample(c(0,1),1)
    X1 <- rep(trt, m*n)
    if(trt==1) { 
      X1[id %% 2 == 0] <- 0 
    }else { 
      X1[id %% 2 == 0] <- 1 
    }
    
    theta <- -(rho*1^2)/(rho-1)
    b_r   <- rep(rnorm(m, 0, theta),n)
    e     <- rnorm(m*n, 0, 1)
    Y     <- beta0 + beta1*X1 + b_r + e
    
    df <- data.frame(id=id, Y=Y, X1=X1)
    
    # GEE
    mod_gee <- gee::gee(Y ~ X1, corstr="exchangeable", id=id, family="gaussian")
    coef_gee <- summary(mod_gee)$coefficients[, c(1,4)]
    inb1_gee <- coef_gee[2,1] - 1.96*coef_gee[2,2] < beta1 & 
                beta1 < coef_gee[2,1] + 1.96*coef_gee[2,2]
    
    # LMM
    mod_lmm <- lme4::lmer(Y ~ X1 + (1|id), REML = TRUE)
    coef_lmm <- summary(mod_lmm)$coefficients[, c(1,2)]
    inb1_lmm <- coef_lmm[2,1] - 1.96*coef_lmm[2,2] < beta1 & 
                beta1 < coef_lmm[2,1] + 1.96*coef_lmm[2,2]
    
    # OLS + FE
    mod_ols <- lm(Y ~ X1 + factor(id))
    coef_ols <- summary(mod_ols)$coefficients[1:2, 1:2]
    inb1_ols <- coef_ols[2,1] - 1.96*coef_ols[2,2] < beta1 & 
                beta1 < coef_ols[2,1] + 1.96*coef_ols[2,2]
    
    # output
    o <- data.frame(model  = c("GEE","LMM","OLS"),
                    b0_est = c(coef_gee[1,1], coef_lmm[1,1], coef_ols[1,1]),
                    b1_est = c(coef_gee[2,1], coef_lmm[2,1], coef_ols[2,1]),
                   #b0_SE  = c(coef_gee[1,2], coef_lmm[1,2], coef_ols[1,2]),
                   #b1_SE  = c(coef_gee[2,2], coef_lmm[2,2], coef_ols[2,2]),
                    inb1   = c(inb1_gee, inb1_lmm, inb1_ols))
    return(o)
}

beta0 <- 1.79
beta1 <- -0.8
res <- tibble(rho=numeric(), m=numeric(), n=numeric(), beta0=numeric(), beta1=numeric(),
              b0_gee=numeric()  , b0_lmm=numeric()  , b0_ols=numeric(),
              b1_gee=numeric()  , b1_lmm=numeric()  , b1_ols=numeric(),
              inb1_gee=numeric(), inb1_lmm=numeric(), inb1_ols=numeric())

for(rho in c(0.1, 0.7)){
  for(m in c(5, 30)) {
    for(n in c(3, 10, 50)) {
      
      sims_p32 <- tibble(model = character(), 
                         b0_est = numeric(), b1_est = numeric(), 
                         #b0_SE = numeric(), b1_SE = numeric(),
                         in_b1 = numeric())
      for(i in 1:200) {
        sims_p32 <- rbind(sims_p32, sim_p32(rho=rho, m=m, n=n, beta0=beta0, beta1=beta1))
      }
      res <- res |> add_row(rho=rho, m=m, n=n, beta0=beta0, beta1=beta1,
                            b0_gee   = beta0 - mean(sims_p32[sims_p32$model=="GEE",]$b0_est),
                            b0_lmm   = beta0 - mean(sims_p32[sims_p32$model=="LMM",]$b0_est),
                            b0_ols   = beta0 - mean(sims_p32[sims_p32$model=="OLS",]$b0_est),
                            b1_gee   = beta1 - mean(sims_p32[sims_p32$model=="GEE",]$b1_est),
                            b1_lmm   = beta1 - mean(sims_p32[sims_p32$model=="LMM",]$b1_est),
                            b1_ols   = beta1 - mean(sims_p32[sims_p32$model=="OLS",]$b1_est),
                            #SEb0_gee = mean(sims_p32[sims_p32$model=="GEE",]$b0_SE),
                            #SEb0_lmm = mean(sims_p32[sims_p32$model=="LMM",]$b0_SE),
                            #SEb0_ols = mean(sims_p32[sims_p32$model=="OLS",]$b0_SE),
                            inb1_gee = mean(sims_p32[sims_p32$model=="GEE",]$inb1),
                            inb1_lmm = mean(sims_p32[sims_p32$model=="LMM",]$inb1),
                            inb1_ols = mean(sims_p32[sims_p32$model=="OLS",]$inb1))      
    }
  }
}
res_p32 <- res
rm(list = c("beta0","beta1","i","m","n","rho","res"))
colnames(res_p32) <- c("rho","m","n","beta0","beta1",   
                       "b0 bias (gee)","b0 bias (lmm)","b0 bias (ols)",
                       "b1 bias (gee)","b1 bias (lmm)","b1 bias (ols)",
                       "b1 covrge (gee)","b1 covrge (lmm)","b1 covrge (ols)")











path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw3/hw3_objects.Rdata"))



################################ D E P R E C A T E D ####################################

### (Create correlated data by hand)
#rm(list = ls())
#sim2 <- function(m=300, n=3, rho=0.25, beta0=-1.5, beta1=0.5, beta2=0.5) {
#  
#  # Create correlated e vector with exchangeable correlation structure
#  Sigma <- toeplitz(c(1,rho,rho))
#  e <- mvtnorm::rmvnorm(n=m, mean=mu, sigma=Sigma) |> t() |> as.vector()
#  id <- rep(1:m, each=n)
#  
#  # Create X1 vector
#  sex <- sample(c(0,1),1)
#  X1 <- rep(sex, m*n)
#  if(sex==1) { 
#    X1[id %% 2 == 0] <- 0 
#  }else { 
#    X1[id %% 2 == 0] <- 1 
#  }
#  
#  # Create X2 vector
#  X2 <- rep(1:n, m)
#  
#  # Create Y vector with true relationship
#  latent_Y <- beta0 + beta1*X1 + beta2*X2 + e 
#  Y <- rbinom(n=length(latent_Y), size=1, prob=invlogit(latent_Y)) #<-try with expit
#  
#  # Pack it all in a dataframe
#  df <- data.frame(id=id, Y=Y, 
#                   X1=matrix(X1, nrow=m*n, ncol=1),
#                   X2 = X2) 
#  
#  # Run GEE
#  mod    <- gee(Y ~ X1 + X2, corstr="exchangeable", id=id, data=df, family="binomial") 
#  modsum <- summary(mod)
#  modrho <- modsum$working.correlation[1,2]
#  mod_b0 <- modsum$coefficients[1,1]
#  mod_se0<- modsum$coefficients[1,4]
#  mod_b1 <- modsum$coefficients[2,1]
#  mod_se1<- modsum$coefficients[2,4]
#  mod_b2 <- modsum$coefficients[3,1]
#  mod_se2<- modsum$coefficients[3,4]
#  
#  o <- list("rho" = modrho, 
#            "b0"  = mod_b0,
#            "se0" = mod_se0,
#            "b1"  = mod_b1,
#            "se1" = mod_se1,
#            "b2"  = mod_b2,
#            "se2" = mod_se2)
#  return(o)
#}

### (First Attempt)
### simul_logit <- function(m, n, rho, beta0, beta1, beta2, sigma=1) {
###  # create correlated X which will determine probabilities p
###  theta = -(rho*sigma^2)/(rho-1)
###  b_r <- rep(rnorm(m, 0, theta),n)
###  
###  # binary gender
###  id <- rep(1:m,n)
###  sex <- sample(c(0,1),1)
###  X1 <- rep(sex, m*n)
###  if(sex==1) { 
###    X1[id %% 2 == 0] <- 0 
###  }else { 
###    X1[id %% 2 == 0] <- 1 
###  }
###  
###  # Time var
###  X2 <- rep(1:n, m)
###  
###  # construct probabilities that determine P(Y=1)
###  p <- exp(beta0+b_r+beta1*X1+beta2*X2)/(1+exp(beta0+b_r+beta1*X1+beta2*X2))
###  Y <- rbinom(n=m*n, size=1, prob=p)
###  id <- rep(1:m, n)
###  df <- data.frame(id=id, Y=Y, 
###                   X1=matrix(X1, nrow=m*n, ncol=1), 
###                   X2=X2,
###                   p=matrix(p, nrow=m*n, ncol=1)) |> dplyr::arrange(id)
###}
