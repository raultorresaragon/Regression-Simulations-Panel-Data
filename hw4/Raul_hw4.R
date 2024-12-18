library(dplyr)
library(tidyr)
library(gee)
library(lme4)
library(SimCorMultRes)
set.seed(571)
rm(list = ls())
df <- readr::read_table("hw4/sixcity.dat", col_names = c("ws","id","age","smokes"))


# ~~~~~~~~~ #
# Problem 1 #
# ~~~~~~~~~ #

# (1)
mod1 <- lme4::glmer(ws ~ age*smokes + (1|id), data=df, family = binomial)
mod1_coef <- summary(mod1)$coefficients


# (2)
mod2 <- gee::gee(ws ~ age*smokes, corstr="exchangeable", id=id, family=binomial, data=df)
mod2_coef <- summary(mod2)$coefficients[,c(1,4,5)]


# (3) BONUS




# ~~~~~~~~~ #
# Problem 2 #
# ~~~~~~~~~ #

# run one sim
# -----------

run_sim <- function(m=30, n=6, rho=0.6, beta0=1, beta1=3, ymethod="rbin") {
  
  # Create correlated e vector with exchangeable correlation structure
  Sigma <- toeplitz(c(1,rep(rho, n-1)))
  e <- mvtnorm::rmvnorm(n=m, mean=rep(0, n), sigma=Sigma) |> t() |> as.vector()
  id <- rep(1:m, each=n)
  
  # Create X1 vector
  trt <- sample(c(0,1),1)
  X1 <- rep(trt, m*n)
  if(trt==1) { 
    X1[id %% 2 == 0] <- 0 
  }else { 
    X1[id %% 2 == 0] <- 1 
  }
  
  # Create Y vector with true relationship
  if(ymethod == "rbin"){
    Y <- rbin(clsize = n, intercepts = beta0, betas = c(beta1), 
              xformula = ~X1, cor.matrix = Sigma, link = "logit")$simdata$y
  } else {
    theta = -(rho*1^2)/(rho-1)
    b_r <- sort(rep(rnorm(m, 0, theta),n))
    p <- exp(beta0+b_r+beta1*X1)/(1+exp(beta0+b_r+beta1*X1))
    Y <- rbinom(n=m*n, size=1, prob=p)
  }
  # Pack it all in a dataframe
  df <- data.frame(id=id, Y=Y, X1=matrix(X1, nrow=m*n, ncol=1))
  
  # fit GEE
  coef_gee <- try(summary(
    gee(Y~X1, corstr="exchangeable", id=id, data=df, family="binomial", silent = TRUE)
  )$coefficients[,c(1,4,5)], silent = TRUE)

  
  # fit GLMM
  coef_glmm <- try(summary(
    lme4::glmer(Y ~ X1 + (1|id), data=df, family = binomial)
  )$coefficients, silent = TRUE) 

  
  if(class(coef_gee)[1]=="matrix" & class(coef_glmm)[1]=="matrix") {
    
    inside_gee <- as.numeric(coef_gee[2,1]-1.96*(coef_gee[2,2]) < beta1 & 
                               beta1 < coef_gee[2,1]+1.96*(coef_gee[2,2]))
    
    inside_glmm <- as.numeric(coef_glmm[2,1]-1.96*(coef_glmm[2,2]) < beta1 & 
                                beta1 < coef_glmm[2,1]+1.96*(coef_glmm[2,2]))
    
    o <- data.frame(m=c(m,m), n=c(n,n), r=c(rho,rho), 
                    beta0=c(beta0,beta0), beta1=c(beta1,beta1),
                    mod=c("gee","glmm"), 
                    b0 = c(coef_gee[1,1], coef_glmm[1,1]),
                    b1 = c(coef_gee[2,1], coef_glmm[2,1]),
                    inside_b1 = c(inside_gee, inside_glmm),
                    se0= c(coef_gee[1,2], coef_glmm[1,2]),
                    se1= c(coef_gee[2,2], coef_glmm[2,2]))
  } else {
    o <- data.frame(m=NA, n=NA, r=NA, beta0=NA, beta1=NA, 
                    mod=NA, b0=NA, b1=NA, inside_b1=NA, se0=NA, se1=NA)
  }
  return(o)
}

# run sim N times
# ---------------

iterate_sim <- function(N=1e2, m=m, n=n, rho=r, beta0=1, beta1=2.5, ymethod="rbin") {
  d <- data.frame(m=numeric(), n=numeric(), r=numeric(), beta0=numeric(), beta1=numeric(), 
                  mod=character(), b0=numeric(), b1=numeric(), 
                  inside_b1 = numeric(), se0=numeric(), se1=numeric())
  for(i in 1:N) {
    d <- rbind(d, run_sim(m=m, n=n, rho=rho, beta0=beta0, beta1=beta1, ymethod=ymethod))
  }
  res <- data.frame(m=m, n=n, r=rho, 
           beta0=c(beta0,beta0), 
           beta1=c(beta1,beta1),
           mod=c("gee","glmm"), 
           b0= c(mean(d[d$mod=="gee",]$b0,na.rm=TRUE), mean(d[d$mod=="glmm",]$b0,na.rm=TRUE)),
           b1= c(mean(d[d$mod=="gee",]$b1,na.rm=TRUE), mean(d[d$mod=="glmm",]$b1,na.rm=TRUE)),
    inside_b1= c(mean(d[d$mod=="gee",]$inside_b1,na.rm=TRUE), mean(d[d$mod=="glmm",]$inside_b1,na.rm=TRUE)),
           se0=c(sd(d[d$mod=="gee",]$b0,na.rm=TRUE), sd(d[d$mod=="glmm",]$b0,na.rm=TRUE)),
           se1=c(sd(d[d$mod=="gee",]$b1,na.rm=TRUE), sd(d[d$mod=="glmm",]$b1,na.rm=TRUE)))
  return(res)
}


# run sims N times per m,n,rho combo
# ----------------------------------

tab3 <- data.frame(m=numeric(), n=numeric(), r=numeric(), 
                   beta0=numeric(), beta1=numeric(), mod=character(), 
                   b0=numeric(), b1=numeric(), se0=numeric(), se1=numeric())

tab3_br <- data.frame(m=numeric(), n=numeric(), r=numeric(), 
                      beta0=numeric(), beta1=numeric(), mod=character(), 
                      b0=numeric(), b1=numeric(), se0=numeric(), se1=numeric())

for(r in c(0.1, 0.8)) {
  for(m in c(20,50,100)) {
    for(n in c(5,10,25)) {
      #tab3 <- rbind(tab3, iterate_sim(N=100, m=m, n=n, rho=r, beta0=1, beta1=2.5))
      tab3_br <- rbind(tab3_br, iterate_sim(N=100, m=m, n=n, rho=r, beta0=1, beta1=2.5, ymethod="br"))
    }
  }
}
colnames(tab3)[colnames(tab3) == "inside_b1"] <- "b1 coverage"
colnames(tab3_br)[colnames(tab3_br) == "inside_b1"] <- "b1 coverage"


# PLOTS
myplot <- function(mydata, myrho, ci=0.9) {
  
  mycols <- c("darkorange","darkgreen")
  z <- qnorm(ci)
  mydata$ci_25 <- mydata$b1-z*mydata$se1
  mydata$ci_975 <- mydata$b1+z*mydata$se1
  offset <- c(-0.03,+0.03,-0.03,+0.03,-0.03,+0.03)
  
  if(myrho == 0.1) {
    myat <- c(1,2,2.5,3,4)
    mylabels <- c("1","2",expression(paste(beta,"_1")),"3","4")
    mytitle <- expression(paste(rho, "=0.1"))
    myxlim <- c(1,4)
  } else {
    myat <- c(0,1,2,2.5,3,4,5)
    mylabels <- c("0","1","2",expression(paste(beta,"_1")),"3","4","5")
    mytitle <- expression(paste(rho, "=0.8"))
    myxlim <- c(0,5)
  }
  
  # TITLE AND LEGEND
  par(mfrow=c(1,1), mar = c(0,0,2,0))
  layout(mat = matrix(c(1,1,1,2,3,4), nrow=2, ncol=3, byrow=TRUE), 
         heights = c(0.2, 0.8))
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="", main = mytitle, cex.main=2.5)
  legend(x = "bottom", inset = 0,
         legend = c("GEE", "GLMM"), 
         col=mycols, 
         lwd=5, cex=2, horiz = TRUE)
  par(mar = c(4,2,2,0.5))
  
  # m=100
  mydata100 <- mydata[mydata$m==100 & mydata$r==myrho, ]
  plot(as.numeric(factor(n))~b1, col=mycols[factor(mod)], data=mydata[mydata$m==100 & mydata$r==myrho, ], 
       xlim=myxlim, 
       xlab=expression(paste(beta, " coefficients")), 
       ylab="n", 
       main="m=100",
       xaxt="n", yaxt="n",
       cex = 1.5)
  axis(side=1, at=myat, labels = mylabels)
  axis(side=2, at=c(1,2,3), labels = c("5","10","25"))
  abline(v=2.5, col="gray", lwd=2, lty=2)
  segments(mydata100$ci_25, as.numeric(factor(mydata100$n)) + offset, 
           mydata100$ci_975, as.numeric(factor(mydata100$n)) + offset, 
           lty=3, lwd=1, col = mycols[factor(mydata100$mod)])
  #points(n~b0, col=mycols[factor(mod)], data=mydata[mydata$m==100 & mydata$r==myrho, ])
  #abline(v=1  , col="gray", lwd=2, lty=2)
  
  # m=50
  mydata50 <- mydata[mydata$m==50 & mydata$r==myrho, ]
  plot(as.numeric(factor(n))~b1, col=mycols[factor(mod)], data=mydata[mydata$m==50 & mydata$r==myrho, ], 
       xlim=myxlim, xlab=expression(paste(beta, " coefficients")), ylab="n", 
       main="m=50", xaxt="n", yaxt="n",
       cex = 1.5)
  axis(side=1, at=myat, labels = mylabels)
  axis(side=2, at=c(1,2,3), labels = c("5","10","25"))
  abline(v=2.5, col="gray", lwd=2, lty=2)
  segments(mydata50$ci_25, as.numeric(factor(mydata100$n)) + offset, 
           mydata50$ci_975, as.numeric(factor(mydata100$n)) + offset, 
           lty=3, lwd=1, col = mycols[factor(mydata50$mod)])  
  #points(n~b0, col=mycols[factor(mod)], data=mydata[mydata$m==50 & mydata$r==myrho, ])
  #abline(v=1  , col="gray", lwd=2, lty=2)
  
  # m=20
  mydata20 <- mydata[mydata$m==20 & mydata$r==myrho, ]
  plot(as.numeric(factor(n))~b1, col=mycols[factor(mod)], data=mydata[mydata$m==20 & mydata$r==myrho, ], 
       xlim=myxlim, xlab=expression(paste(beta, " coefficients")), ylab="n", 
       main="m=20", xaxt="n", yaxt="n",
       cex = 1.5)
  axis(side=1, at=myat, labels = mylabels)
  axis(side=2, at=c(1,2,3), labels = c("5","10","25"))
  abline(v=2.5, col="gray", lwd=2, lty=2)
  segments(mydata20$ci_25, as.numeric(factor(mydata100$n)) + offset, 
           mydata20$ci_975, as.numeric(factor(mydata100$n)) + offset, 
           lty=3, lwd=1, col = mycols[factor(mydata20$mod)])  
  #points(n~b0, col=mycols[factor(mod)], data=mydata[mydata$m==20 & mydata$r==myrho, ])
  #abline(v=1  , col="gray", lwd=2, lty=2)
}

myplot(tab3, 0.1)
myplot(tab3, 0.8)
myplot(tab3_br, 0.1)
myplot(tab3_br, 0.8)

path <- "/Users/raulta/Library/CloudStorage/GoogleDrive-rauldavidt@gmail.com/My Drive"
save.image(file = paste0(path,"/_STAT 571 Advanced Regression/hw4/hw4_objects.Rdata"))

