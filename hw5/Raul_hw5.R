library(lme4)
library(SimCorMultRes)
library(gee)
library(dplyr)
set.seed(571)
rm(list = ls())

# Problem 1
# m = 25
# n = 6
# rho = 0.61
# beta0 = 1
# beta1 = 2.5
# beta2 = 0.5
# miss_freq = 0.3

# ~~~~~~~~~~~~~~~~~~~~~~ #
# Create correlated data #
# ~~~~~~~~~~~~~~~~~~~~~~ #

make_corr_df <- function(m, n, rho, beta0, beta1, beta2) {
  id <- rep(1:m, each=n)
  
  # Create X vectors
  trt <- sample(c(0,1),1)
  X1 <- rep(trt, m*n)
  if(trt==1) { 
    X1[id %% 2 == 0] <- 0 
  }else { 
    X1[id %% 2 == 0] <- 1 
  }
  X2 <- rnorm(m*n, 5, 1)
  
  # Create Y vector with true relationship (LMM)
  theta <- -(rho*1^2)/(rho-1)
  b_r <- sort(rep(rnorm(m, 0, theta),n))
  e <- rnorm(n=m*n,0,1)
  Y1 = beta0+b_r+beta1*X1+beta2*X2+e
  
  # Create Y vector with true relationship (GEE exchangeability)
  Sigma <- toeplitz(c(1,rep(rho, n-1)))
  e <- mvtnorm::rmvnorm(n=m, mean=rep(0, n), sigma=Sigma) |> t() |> as.vector()
  Y2 = beta0+beta1*X1+beta2*X2+e
  
  
  # Pack it all in a dataframe
  df <- data.frame(id=id, Y1=Y1, Y2=Y2, X1=matrix(X1, nrow=m*n, ncol=1), X2=X2) |> 
        group_by(id) |> mutate(n = row_number())
  df
}

# ~~~~~~~~~~~~~~~~~~~~~ #
# Introduce missingness #
# ~~~~~~~~~~~~~~~~~~~~~ #

make_NA <- function(index, df) {
  df$fakeY <- NA
  df$fakeY[index] <- -99 #CHANGE TO IF ELSE??
  df$fakeY[df$n == 1] <- NA
  df <- df |> group_by(id) |> tidyr::fill(fakeY, .direction = "down")
  df$Y[df$fakeY == -99] <- NA
  df <- df |> dplyr::select(-n, -fakeY)
}

make_miss_df <- function(df, miss_freq = 0.3, m, n, Y="Y1") {
  
  colnames(df)[colnames(df)==Y] <- "Y"
  
  # MCAR
  Nt <- m*n
  p_miss <- miss_freq
  index <- which(rbinom(n=m*n, size=1, prob = p_miss) == 1)
  df_mcar <- make_NA(index, df)
  
  # MAR 'males' more likely to be missing
  p_miss <- if_else(df$X1==1, miss_freq*2, 0.04)
  index <- which(rbinom(n=length(p_miss), size=1, prob = p_miss) == 1)
  df_mar <- make_NA(index, df)
  
  # MNAR "high Y= very likely missing"
  p_miss <- if_else(df$Y < quantile(df$Y, miss_freq), 0.9, 0.04)
  index  <- which(rbinom(n=length(p_miss), size=1, prob = p_miss) == 1)
  df_mnar <- make_NA(index = index, df)
  
  o <- list(df_mcar = df_mcar, df_mar = df_mar, df_mnar = df_mnar)
  o
  
}

# ~~~~~~~~~~ #
# FIT MODELS #
# ~~~~~~~~~~ #

fit_mods <- function(df, m, n, rho, beta1) {
  
  df <- df[complete.cases(df), ]
  
  # fit GEE exchangeable
  coef_gee <- try(summary(
    gee(Y~X1+X2, corstr="exchangeable", id=id, data=df, silent = TRUE)
  )$coefficients[,c(1,4,5)], silent = TRUE)
  
  
  # fit GEE AR(1)
  coef_gee_un <- try(summary(
    gee(Y~X1+X2, corstr="independence", id=id, data=df, silent = TRUE)
  )$coefficients[,c(1,4,5)], silent = TRUE)
  
  
  # fit LMM
  coef_lmm <- try(summary(
    lme4::lmer(Y ~ X1 + X2 + (1|id), data=df)
  )$coefficients, silent = TRUE) 
  
  
  if(class(coef_gee)[1]=="matrix" & 
     class(coef_lmm)[1]=="matrix" & 
     class(coef_gee_un)[1]=="matrix") {
    
    inside_gee <- as.numeric(coef_gee[2,1]-1.96*(coef_gee[2,2]) < beta1 & 
                               beta1 < coef_gee[2,1]+1.96*(coef_gee[2,2]))
    
    inside_gee_un<-as.numeric(coef_gee_un[2,1]-1.96*(coef_gee_un[2,2]) < beta1 & 
                               beta1 < coef_gee_un[2,1]+1.96*(coef_gee_un[2,2]))
    
    inside_lmm <- as.numeric(coef_lmm[2,1]-1.96*(coef_lmm[2,2]) < beta1 & 
                               beta1 < coef_lmm[2,1]+1.96*(coef_lmm[2,2]))
    
    signif_gee <- as.numeric(!(coef_gee[2,1]-1.96*(coef_gee[2,2]) < 0 & 
                             0 < coef_gee[2,1]+1.96*(coef_gee[2,2])))
    
    signif_gee_un<-as.numeric(!(coef_gee_un[2,1]-1.96*(coef_gee_un[2,2]) < 0 & 
                             0 < coef_gee_un[2,1]+1.96*(coef_gee_un[2,2])))
    
    signif_lmm <- as.numeric(!(coef_lmm[2,1]-1.96*(coef_lmm[2,2]) < 0 & 
                             0 < coef_lmm[2,1]+1.96*(coef_lmm[2,2])))    
    
    o <- data.frame(m=m, n=n, r=rho, 
                    pct_miss= 1-nrow(df)/(m*n),
                    n_cases=length(unique(df$id)),
                    beta1=beta1, 
                    mod=c("lmm","gee_EX","gee_UN"), 
                    b1 = c(coef_lmm[2,1], coef_gee[2,1], coef_gee_un[2,1]),
                    inside_b1 = c(inside_lmm, inside_gee, inside_gee_un),
                    signif_b1 = c(signif_lmm, signif_gee, signif_gee_un),
                    se1= c(coef_lmm[2,2], coef_gee[2,2], coef_gee_un[2,2]))
                    
  } else {
    o <- data.frame(m=NA, n=NA, r=NA, pct_miss=NA, n_cases=NA, 
                    beta1=NA, mod=NA, b1=NA, inside_b1=NA, signif_b1=NA, se1=NA)
  }
  o
}


# ~~~~~~~ #
# RUN SIM #
# ~~~~~~~ #

run_sim <- function(m, n, rho, beta0, beta1, beta2, miss_freq, Y_type = "Y1") {
  
  # create correlated df
  df <- make_corr_df(m, n, rho, beta0, beta1, beta2) 
  
  # create missing dfs
  miss_dfs <- make_miss_df(df, miss_freq, m, n, Y_type)
  df_mcar <- miss_dfs$df_mcar
  df_mar  <- miss_dfs$df_mar
  df_mnar <- miss_dfs$df_mnar
  
  # fit models on available cases
  fit_mcar <- fit_mods(df_mcar, m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mcar")
  fit_mar  <- fit_mods(df_mar,  m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mar")
  fit_mnar <- fit_mods(df_mnar, m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mnar")
  stopifnot(ncol(fit_mcar) == ncol(fit_mar) & ncol(fit_mnar) == ncol(fit_mar))
  comp_cases <- rbind(fit_mcar, fit_mar, fit_mnar)
  comp_cases$cases <- "available"
  
  # fit models on compelte cases: drop clusters with missingness
  df_mcar <- df_mcar |> group_by(id) |> mutate(has_NA = max(is.na(Y))) |> filter(has_NA==0)
  df_mar  <- df_mar  |> group_by(id) |> mutate(has_NA = max(is.na(Y))) |> filter(has_NA==0)
  df_mnar <- df_mnar |> group_by(id) |> mutate(has_NA = max(is.na(Y))) |> filter(has_NA==0)
  fit_mcar <- fit_mods(df_mcar, m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mcar")
  fit_mar  <- fit_mods(df_mar,  m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mar")
  fit_mnar <- fit_mods(df_mnar, m=m, n=n, r=rho, beta1=beta1) |> mutate(miss_type = "mnar")
  stopifnot(ncol(fit_mcar) == ncol(fit_mar) & ncol(fit_mnar) == ncol(fit_mar))
  avai_cases <- rbind(fit_mcar, fit_mar, fit_mnar)
  avai_cases$cases <- "complete"
  
  #stopifnot(ncol(comp_cases) == ncol(avai_cases))
  o <- rbind(comp_cases, avai_cases)
  return(o)
  
}

iterate_sim <- function(N=1e2, m, n, rho, beta0, beta1, beta2, miss_freq) {
  d <- data.frame(m=numeric(), n=numeric(), r=numeric(), n_cases=numeric(),
                  pct_miss=numeric(), 
                  beta1=numeric(), mod=character(), b1=numeric(), 
                  inside_b1 = numeric(), 
                  signif_b1 = numeric(), 
                  se1=numeric(),
                  miss_type = character(), 
                  cases = character())
  for(i in 1:N) {
    d <-rbind(d,run_sim(m, n, rho, beta0, beta1, beta2, miss_freq=miss_freq, Y_type="Y1"))
  }
  res <- d |> group_by(m, n, r, mod, miss_type, cases) |> 
    reframe(pct_miss = mean(pct_miss, na.rm=TRUE),
            n_cases = max(n_cases),
            beta1 = mean(beta1), 
            b1_bias = beta1 - mean(b1, na.rm=TRUE),
            coverage = mean(inside_b1),
            signif_rate = mean(signif_b1, na.rm=TRUE),
            emp_se1 = sd(se1, na.rm=TRUE))
                             
  res <- res |> 
         dplyr::select(m, n, r, beta1, pct_miss, miss_type, 
                       n_cases, cases, mod, b1_bias, coverage, signif_rate, emp_se1)
  return(res)
}

tab1 <- iterate_sim(N=100, m=50, n=6, rho=0.68, 
                    beta0=1, beta1=2, beta2=0.5, miss_freq=0.3) |> 
        dplyr::arrange(m,n,r,beta1, cases, miss_type, mod)

tab2 <- iterate_sim(N=100, m=30, n=10, rho=0.68, 
                    beta0=1, beta1=2, beta2=0.5, miss_freq=0.3) |> 
        dplyr::arrange(m,n,r,beta1, cases, miss_type, mod)

tab3 <- iterate_sim(N=100, m=65, n=4, rho=0.68, 
                    beta0=1, beta1=2, beta2=0.5, miss_freq=0.3) |> 
        dplyr::arrange(m,n,r,beta1, cases, miss_type, mod)


beautify_tab <- function(tab) {
  tab <- tab |> dplyr::filter(!is.na(m)) |> 
                tidyr::pivot_wider(id_cols = c("m","n","r","beta1","pct_miss",
                                               "miss_type","n_cases","cases"),
                                   values_from = c("b1_bias", "coverage","signif_rate"),
                                   names_from = c("mod"))
  colnames(tab) <- c("m","n","r", "beta1", "prop missing", 
                     "miss type", "n cases", "cases",
                     "bias (GEE EX)", "bias (GEE UN)", "bias (LMM)",
                     "coverage (GEE EX)", "coverage (GEE UN)", "coverage (LMM)",
                     "signif rate (GEE EX)", "signif rate (GEE UN)", "signif rate (LMM)")
  return(tab)
}

tab1 <- beautify_tab(tab1)
tab2 <- beautify_tab(tab2)
tab3 <- beautify_tab(tab3)

make_bigtab <- function(N=100, m=26, mf = 0.3) {
  bigtab <- rbind(
         iterate_sim(N, m=m, n= 5, rho=0.61, beta0=1, beta1=2, beta2=0.5, miss_freq=mf),
         iterate_sim(N, m=m, n=10, rho=0.61, beta0=1, beta1=2, beta2=0.5, miss_freq=mf),
         iterate_sim(N, m=m, n=20, rho=0.61, beta0=1, beta1=2, beta2=0.5, miss_freq=mf),
         iterate_sim(N, m=m, n=40, rho=0.61, beta0=1, beta1=2, beta2=0.5, miss_freq=mf),
         iterate_sim(N, m=m, n=80, rho=0.61, beta0=1, beta1=2, beta2=0.5, miss_freq=mf)
       ) |> filter(!is.na(m)) |> 
            mutate(prop_missing = as.character(mf+0.1)) |> 
            dplyr::select(-pct_miss)
}
hugetab <- rbind(make_bigtab(mf=0.1), make_bigtab(mf=0.2), make_bigtab(mf=0.3)) 


mycols <- c("darkorange","darkgreen","darkred")
make_subplot <- function(mi_type, hugetab, mf) {
 plotdata <- hugetab[hugetab$cases=="available" & 
                      hugetab$miss_type==mi_type & 
                      hugetab$prop_missing==mf, ] |> 
             arrange(mod, n)
 plot(b1_bias~n, col=mycols[factor(mod)], data=plotdata[plotdata$miss_type == mi_type, ],
     main = stringr::str_to_upper(mi_type),
     xlab = "n",
     ylab = paste0("proportion missing = ", mf),
     xaxt = "n",
     ylim = c(-0.55, 0.65))
 axis(side=1, at = c(5,10,20,40,60,80), labels = c(5,10,20,40,60,80))
 lines(b1_bias~n, col="darkorange", data=plotdata[plotdata$mod == "gee_EX", ])
 lines(b1_bias~n, col="darkgreen",data=plotdata[plotdata$mod == "gee_UN", ])
 lines(b1_bias~n, col="darkred",data=plotdata[plotdata$mod == "lmm", ])
 abline(h = 0, lty = 3)
}

make_plot <- function(mf) {

  # TITLE AND LEGEND
  par(mfrow=c(1,1), mar = c(0,0,2,0))
  layout(mat = matrix(c(1,1,1,2,3,4), nrow=2, ncol=3, byrow=TRUE), 
         heights = c(0.2, 0.8))
  #layout(mat = matrix(c(1,1,1,2,3,4,5,6,7,8,9,10), nrow=4, ncol=3, byrow=TRUE), 
  #       heights = c(0.3, 0.7, 0.7, 0.7))  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="", 
       main = bquote(paste(hat(beta)[1], " bias with ", .(100*mf), "% missingness")), 
       cex.main=1.8)
  legend(x = "bottom", inset = 0,
         legend = c("GEE (Ex)", "GEE (Un)", "LMM"), 
         col=mycols, 
         lwd=5, cex=1.5, horiz = TRUE)
  
  # SUBPLOTS
  par(mar = c(5,2,2,0.5))
  make_subplot("mar",  hugetab, mf = mf)
  make_subplot("mcar", hugetab, mf = mf)
  make_subplot("mnar", hugetab, mf = mf)
  
  #make_subplot("mar",  hugetab, mf = 0.3)
  #make_subplot("mcar", hugetab, mf = 0.3)
  #make_subplot("mnar", hugetab, mf = 0.3)
  #
  #make_subplot("mar",  hugetab, mf = 0.4)
  #make_subplot("mcar", hugetab, mf = 0.4)
  #make_subplot("mnar", hugetab, mf = 0.4)

}
make_plot(mf=0.1)
make_plot(mf=0.2)
make_plot(mf=0.3)
make_plot(mf=0.4)



# ~~~~~~~~~ #
# Problem 3 #
# ~~~~~~~~~ #

sim_p3 <- function(m, n, beta0, beta1, lambda = 1) {
  
  # create cluster id
  id <- rep(1:m, each=n)
  
  # Create X vectors
  trt <- sample(c(0,1),1)
  X1 <- rep(trt, m*n)
  if(trt==1) { 
    X1[id %% 2 == 0] <- 0 
  }else { 
    X1[id %% 2 == 0] <- 1 
  }

  # Create Y vector with true relationship (LMM)
  # theta <- -(rho*1^2)/(rho-1)
  b_r <- sort(rep(rexp(n=m, lambda),n))
  e <- rnorm(n=m*n,0,2)
  Y = beta0+b_r+beta1*X1+e

  # Pack it all in a dataframe
  df <- data.frame(id=id, Y=Y, X1=matrix(X1, nrow=m*n, ncol=1)) |> 
        group_by(id) |> mutate(n = row_number())

  # fit GEE exchangeable
  coef_gee <- try(summary(
    gee(Y~X1, corstr="exchangeable", id=id, data=df, silent = TRUE)
  )$coefficients[,c(1,4,5)], silent = TRUE)
  
  # fit LMM
  coef_lmm <- try(summary(
    lme4::lmer(Y ~ X1 + (1|id), data=df)
  )$coefficients, silent = TRUE) 


  if(class(coef_gee)[1]=="matrix" & class(coef_lmm)[1]=="matrix") {
    
    inside_gee <- c(as.numeric(coef_gee[1,1]-1.96*(coef_gee[1,2]) < beta0 & 
                               beta0 < coef_gee[1,1]+1.96*(coef_gee[1,2])),
                    as.numeric(coef_gee[2,1]-1.96*(coef_gee[2,2]) < beta1 & 
                                 beta1 < coef_gee[2,1]+1.96*(coef_gee[2,2])))
    
    inside_lmm <- c(as.numeric(coef_lmm[1,1]-1.96*(coef_lmm[1,2]) < beta0 & 
                               beta0 < coef_lmm[1,1]+1.96*(coef_lmm[1,2])),
                    as.numeric(coef_lmm[2,1]-1.96*(coef_lmm[2,2]) < beta1 & 
                                 beta1 < coef_lmm[2,1]+1.96*(coef_lmm[2,2])))
    
    signif_gee <- as.numeric(!(coef_gee[2,1]-1.96*(coef_gee[2,2]) < 0 & 
                                 0 < coef_gee[2,1]+1.96*(coef_gee[2,2])))
    
    signif_lmm <- as.numeric(!(coef_lmm[2,1]-1.96*(coef_lmm[2,2]) < 0 & 
                                 0 < coef_lmm[2,1]+1.96*(coef_lmm[2,2])))    
    
    o <- data.frame(m=m, n=n, beta1=beta1, l=lambda, 
                    mod       = c("lmm", "gee_EX"), 
                    b0        = c(coef_lmm[1,1], coef_gee[1,1]),
                    b1        = c(coef_lmm[2,1], coef_gee[2,1]),
                    inside_b1 = c(inside_lmm[2], inside_gee[2]),
                    inside_b0 = c(inside_lmm[1], inside_gee[1]),
                    signif_b1 = c(signif_lmm, signif_gee),
                    se1       = c(coef_lmm[2,2], coef_gee[2,2]))
    
  } else {
    o <- data.frame(m=NA, n=NA, beta1=NA, l=NA,
                    mod=NA, b0=NA, b1=NA, 
                    inside_b0=NA, inside_b1=NA, signif_b1=NA, se1=NA)
  }
  o
}

iterate_sim3 <- function(N=1e2, m, n, beta0, beta1, lambda) {
  d <- data.frame(m=numeric(), n=numeric(), l=numeric(), beta1=numeric(), mod=character(), 
                  b0=numeric(),
                  b1=numeric(), 
                  inside_b0 = numeric(),
                  inside_b1 = numeric(), 
                  signif_b1 = numeric(), 
                  se1=numeric())
  for(i in 1:N) {
    d <-rbind(d, sim_p3(m, n, beta0, beta1, lambda))
  }
  res <- d |> group_by(m, n, l, mod) |> 
    reframe(beta1 = mean(beta1), 
            b0_bias = beta0 - mean(b0, na.rm = TRUE),
            b1_bias = beta1 - mean(b1, na.rm=TRUE),
            coverage_b0 = mean(inside_b0),
            coverage_b1 = mean(inside_b1),
            signif_b1 = mean(signif_b1, na.rm=TRUE),
            emp_se1 = sd(se1, na.rm=TRUE))
  
  res <- res |> 
    dplyr::select(m, n, l, beta1, mod, 
                  b0_bias, b1_bias, coverage_b0, coverage_b1, signif_b1, emp_se1)
  return(res)
}

d3 <- data.frame(m=numeric(), n=numeric(), l=numeric(), beta1=numeric(), mod=character(), 
                b0=numeric(), b1=numeric(), 
                inside_b0=numeric(), inside_b1=numeric(), 
                signif_b1=numeric(), 
                se1=numeric())
for(l in c(0.5, 1, 2.5)) {
  for(m in c(5,20,100)) {
    for(n in c(5,20,100)){
      if(m != n) {
        cat(paste0("m=",m, " n=", n, " l=",l))
        d3 <- rbind(d3, iterate_sim3(N=100, m=m, n=n, beta0=1, beta1=2.4, lambda=l))
      }
    }
  }
}

View(d3)
      
path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw5/hw5_objects.Rdata"))



####  D E P R E C A T E D  ####

## d <- data.frame(m=numeric(), n=numeric(), r=numeric(), beta1=numeric(), 
##                 pct_miss=numeric(), miss_type=character(), 
##                 n_cases=numeric(), cases=character(), mod=character(), 
##                 b1_bias=numeric(), coverage=numeric(), emp_se1=numeric())

### MAR
###glm(X1~X2, data=df, family="binomial") |> predict(type="response") |> unname()

### MNAR
###sort_Y1 <- sort(df$Y1, decreasing = TRUE)
###nmar1 <- sort_Y1[ceiling(miss_freq*m*n)]
###index1 <- which(df$Y1 > nmar1)


## Create Y vector with true relationship
## if(ybinary == FALSE) {
##   theta <- -(rho*sigma^2)/(rho-1)
##   b_r <- rnorm(0, theta)
##   Y = beta0+b_r+beta1*X1+beta2*X2+rnorm(0,1)
## } else{ 
##   if(ymethod == "rbin"){
##     Y <- rbin(clsize = n, intercepts = beta0, betas = c(beta1, beta2), 
##               xformula = ~X1+X2, cor.matrix = Sigma, link = "logit")$simdata$y
##   } else {
##     theta = -(rho*1^2)/(rho-1)
##     b_r <- sort(rep(rnorm(m, 0, theta),n))
##     p <- exp(beta0+b_r+beta1*X1+beta2*X2)/(1+exp(beta0+b_r+beta1*X1+beta2*X2))
##     Y <- rbinom(n=m*n, size=1, prob=p)
##   }
## }
## data.frame(m=m, n=n, r=rho, 
## beta1=beta1,
## mod=c("gee","glmm"), 
## b0= c(mean(d[d$mod=="gee",]$b0,na.rm=TRUE), mean(d[d$mod=="glmm",]$b0,na.rm=TRUE)),
## b1= c(mean(d[d$mod=="gee",]$b1,na.rm=TRUE), mean(d[d$mod=="glmm",]$b1,na.rm=TRUE)),
## inside_b1= c(mean(d[d$mod=="gee",]$inside_b1,na.rm=TRUE), mean(d[d$mod=="glmm",]$inside_b1,na.rm=TRUE)),
## se0=c(sd(d[d$mod=="gee",]$b0,na.rm=TRUE), sd(d[d$mod=="glmm",]$b0,na.rm=TRUE)),
## se1=c(sd(d[d$mod=="gee",]$b1,na.rm=TRUE), sd(d[d$mod=="glmm",]$b1,na.rm=TRUE)))



## (1) 

## Who are the Taylor Swift fans? Are they Millennials, Gen Z, people with high GPA, people in committed relationships?
##   
##   We want to identify who these Taylor Swift fans are. A lot of individuals claim to be Taylor Swift fans, but in reality they state that get more dates or to hide their taste in music. In order to separate the phonies from the real fans, we will administer the Taylor-Swift fandom questionnaire (TSQ)--a carefully written survey developed by the best minds of the 21st century in social science, psychology (marketing and MBA individuals had no say in the questions of this survey) to 200 self-proclaimed participants. Anybody with a consistent score (over two years) higher then 70 out of 100 will be deemed a true fan.  
## 
## Statistical Plan
## 
## We will use a linear mixed model to measure the level of association between age, region of origin (international, West, South, Northeast), binarized gender, relationship status, university (UC Berkeley vs UT Austin), time, academic performance in the form of GPA, self-reported anxiety score and TSQ score over time among US college students.  
## 
## We will interact time with each demographic because we think that fandom may vary overtime differently for different students based on origin, university, gender, and relationship status.  
## 
## Note that others have studied the association between Taylor Swift fandom on academics and mental health. (Lazy McLazyperson etal. JASA 2021). However, their study did not account for self-selection given that all the participants volunteered after hearing about the study at a recent concert he attended.  
## 
## Assumptions
## 
## Because the Taylor-Swift survey is designed with self-proclaimed fans in mind (i.e. only those who truly are obsessed with Taylor Swift will get all questions right), we will assume scores will be normally distributed. We assume that anxiety-free single freshen women from the West with moderate to high GPA will score highest. Hence this will be our baseline profile. We will assume a random intercept where these random effects come from a normal distribution with mean zero.   
## 
## Power assumptions  
## 
## Because we only have 200 participants from these two universities, we assume that we will have enough variation in this sample to detect the above-mentioned controls and Taylor-Swift fandom. Thanks to the wide diversity at these two large universities (UC Berkeley and UT Austin) we believe the West, South, Midwest, and Northeast will be represented to detect the level of association. Assuming, for example, that students from the Midwest vary widely in their level of Taylor-Swift fandom from students from the West (Cohen d=0.45), we should be able to detect a difference with 80% Power and $\alpha = 0.10$.  
## 
## Because people move, get deported, or become Katy Perry fans, there is a potential risk of attrition. Even after an assumed 10% attrition rate, the effects are strong enough to still detect differences across different groups with 180 people.  
## 
## These may seem strong power and sample size assumptions, but a recent study by Daniel Paydarfar etal (JSM upcoming presentation) shows significant effects ($\alpha = 0.05$) with 167 Beyonce fans.  




