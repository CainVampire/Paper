# clear the existing environment varable.
rm(list=ls())

# generate data -----------------------------------------------------------
require(fGarch)

tau <- 0.05 ##[to do] tau=0.01
q1 <- qsstd(tau, mean = 0, sd = 2, nu = 6.35, xi = 0.92)
beta0 <- c(-0.30, 0.9, -0.033, -0.066)

#link function
lx <- function(x){
  (0.09-0.02*x*(x<0)+ 0.01*x*(x>=0))*q1
}

simula=function(n) {
  xsamp1=matrix(0,n,300)
  xsamp2=matrix(0,n,300)
  for (j in 1:300) {
    set.seed(j)
    eps <- rsstd(n, mean = 0, sd = 2, nu = 6.35, xi = 0.92)
    r=numeric(n)
    r[1] = eps[1]
    q <- numeric(n)
    q[1] <- q1
    for(t in 2:n){
      r[t] <- ((0.9*q[t-1] + lx(r[t-1]))/ q1) * eps[t]
      q[t] <- 0.9 * q[t-1] + lx(r[t-1])
    }
    xsamp1[,j]=r
    xsamp2[,j]=q
  }
  return(list(qdat=xsamp2,rdat=xsamp1))
}

n <- 250 #[to do] n=500,1000
data1 <- simula(n)

# parameter estimating ----------------------------------------------------

# Linear Regression Spline with single knot at 0
library(splines)

LRS_est <- function(qdatj, rdatj, j){
  q_ex <- qdatj[,j] 
  r_ex <- rdatj[,j]
  Rmat <- cbind(1, pmax(0, r_ex), pmax(0,-r_ex))
  mx <- ncol(Rmat)
  set.seed(j)
  theta <- 0.9 + runif(1, -0.05, 0.05)
  n <- 250 #[to do] n=500,1000
  Xfun <- function(j=1){
    Xj <- numeric(n)
    for(t in 2:n){ 
      for(i in 1:(t-1)){
        Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i, j]
      }
    }
    return(Xj)
  }
  
  X <- sapply(1:mx, Xfun) 
  y <- numeric(n)
  for(t in 2:n){
    y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
  }
  
  library(quantreg)
  dat <- data.frame(y=y[-1], x=X[-1,])
  names(dat) <- c('y', 'x1', 'x2', 'x3')
  rq1 <- rq(y~.+0, tau=0.05, data=dat)
  yhat <- rq1$fitted.values
  qhat <- yhat + theta^(1:(n-1))*q_ex[1]
  AAD_LRS <- mean(abs(qhat - q_ex[2:n]))
  betahat_LRS <- c(rq1$coefficients[1],theta,rq1$coefficients[-1])
  names(betahat_LRS) <- paste0('beta',1:4)
  return(c(betahat_LRS, 'AAD_LRS'=AAD_LRS))
}

res1 <- sapply(1:300, function(j) LRS_est(data1$qdat, data1$rdat, j))
apply(res1, 1, sd)
apply(res1, 1, mean) - c(beta0, 0)


# CAViaR of Asymmetric Slop Model
CAViAR_est <- function(qdatj, rdatj, j){
  q_ex <- qdatj[,j] 
  r_ex <- rdatj[,j]
  rt1 <- pmax(0, r_ex)
  rt2 <- pmax(0, -r_ex)
  dat2 <- data.frame(qt=q_ex[-1], qt1 = q_ex[-n], rt1=rt1[-n], rt2=rt2[-n])
  
  rq2 <- rq(qt~., tau=0.05, data=dat2)#[to do]tau=0.01
  qhat_CAViaR <- rq2$fitted.values
  AAD_CAViaR <- mean(abs(qhat_CAViaR - dat2$qt))
  betahat_CAViaR <- rq2$coefficients
  names(betahat_CAViaR) <- paste0('beta',1:4)
  bias_LRS <- betahat_CAViaR - beta0
  return(c(betahat_CAViaR, 'AAD_CAViaR'=AAD_CAViaR))
}

res2 <- sapply(1:300, function(j) LRS_est(data1$qdat, data1$rdat, j))
apply(res2, 1, mean) - c(beta0, 0)
apply(res2, 1, sd)


# selection # of knots ------------------------------------------------------------
library(splines)
n <- 1500
data2 <- simula(n) ##[to do] 第2种连接函数l(x)

#check function
rhoFun <- function(u, tau){
  tau*u*(u>0) - (1-tau)*u*(u<0)
}

#GC,SIC,AIC
GC <- function(ytest, yhatest, tau){
  mean(rhoFun(ytest-yhatest, tau))
}

SIC <- function(y, yhat, tau, nknots, degree){
  N <- length(y)
  R <- nknots + degree + 1
  log(mean(rhoFun(y-yhat, tau))) + 1/2*log(N)/N*R
}

AIC <- function(y, yhat, tau, nknots, degree){
  N <- length(y)
  R <- nknots + degree + 1
  log(mean(rhoFun(y-yhat, tau))) + 2*(R+1)/(N-(R+2))
}

knots_list <- list(c(0),
                   c(0,3),
                   c(-10,0,3),
                   c(-10,0,3,5),
                   c(-6,-3,0,3,6))
#attr(bs(r_ex,df=4, degree = 1), 'knots')

select_knots_Fun <- function(knots, tau=0.05, degree=1, method='GC'){#[to do] tau=0.01
  nknots = length(knots)
  Rmat <- bs(r_ex, knots=knots, degree = degree, intercept = T)
  #Rmat <- ns(r_ex, knots=knots, intercept = T) #for Natural Cubic Spline
  mx <- ncol(Rmat)
  n <- nrow(Rmat)
  theta <- 0.9
  cat('evaluation is beginging---------------\n')
  Xfun2 <- function(j=1){
    Xj <- numeric(n)
    for(t in 2:n){
      for(i in 1:(t-1)){
        Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i, j]
      }
    }
    return(Xj)
  }
  X <- sapply(1:mx, Xfun2)
  y <- numeric(n)
  for(t in 2:n){
    y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
  }
  
  library(quantreg)
  ntr <- 1000
  dat <- data.frame(y=y[2:ntr], x=X[2:ntr,])
  names(dat) <- c('y', paste0('x', 1:mx))
  rq1 <- rq(y~.+0, tau= tau, data=dat)
  yhat <- rq1$fitted.values
  ytest <- y[(ntr+1):n]
  X_test_dat <- data.frame(X[(ntr+1):n,])
  names(X_test_dat) <- paste0('x', 1:mx)
  yhatest <- predict(rq1, X_test_dat)
  qhat <- yhat + theta^(1:(ntr-1))*q_ex[1]
  AAD_LRS <- mean(abs(qhat - q_ex[2:ntr]))
  
  CIV <- switch(method, 
                GC = GC(ytest, yhatest, tau),
                SIC = SIC(y[2:ntr], yhat, tau, nknots, degree),
                AIC = AIC(y[2:ntr], yhat, tau, nknots, degree))
  cat('evaluation is finished!------------------------- \n')
  return(c('AAD'=AAD_LRS, 'CIV'=CIV))
}

knots_eval <- function(method='GC', tau=0.05){#[to do] tau=0.01
  AC <- sapply(knots_list, function(knots) select_knots_Fun(knots,
                                        tau=tau,degree=1,method=method))
  jknots <- which.min(AC[2,])
  AAD <- AC[1, jknots]
  return(c('AAD'=AAD, 'nKnots'=jknots))
}

N <- 300
AN <- matrix(nrow=N, ncol=2)
for(group in 1:N){
  q_ex <- data2$qdat[,group]
  r_ex <- data2$rdat[,group]
  AN[group,] <- knots_eval() #[to do]method需要在knots_eval中改变默认值
}
colnames(AN) <- c('AAD', 'Nknots')
apply(AN, 2, mean)
sd(AN[,'AAD'])


# Selection of lambda -----------------------------------------------------

##GC,SIC,AIC
GC1 <- function(ytest, yhatest, tau){
  mean(rhoFun(ytest-yhatest, tau))
}

SIC1 <- function(y, yhat, tau, R){
  N <- length(y)
  log(mean(rhoFun(y-yhat, tau))) + 1/2*log(N)/N*R
}

AIC1 <- function(y, yhat, tau, R){
  N <- length(y)
  log(mean(rhoFun(y-yhat, tau))) + 2*(R+1)/(N-(R+2))
}

select_lambda_Fun <- function(lambda, tau=0.05,method='GC'){
  theta <- 0.9
  n <- length(r_ex)
  cat('evaluation is beginging---------------\n')
  y <- numeric(n)
  for(t in 2:n){
    y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
  }
  
  library(quantreg)
  ntr <- 1000
  dat <- data.frame(y=y[2:ntr], x=r_ex[1:(ntr-1)])
  rq1 <- rqss(y~ qss(x, constraint = 'N', lambda = lambda),tau= tau, data=dat)
  suma <- summary(rq1)
  yhat <- predict(rq1, dat)
  X_test_dat <- data.frame(x=r_ex[(ntr):(n-1)])
  ytest <- y[(ntr+1):n]
  yhatest <- predict(rq1, X_test_dat)
  qhat <- yhat + theta^(1:(ntr-1))*q_ex[1]
  AAD_LRS <- mean(abs(qhat - q_ex[2:ntr]))
  p_lambda <- suma$edf
  CIV <- switch(method, 
                GC = GC1(ytest, yhatest, tau),
                SIC = SIC1(y[2:ntr], yhat, tau, p_lambda),
                AIC = AIC1(y[2:ntr], yhat, tau, p_lambda))
  cat('evaluation is finished!------------------------- \n')
  p_lambda
  lambda
  return(c('AAD'=AAD_LRS, 'CIV'=CIV, 'lambda'=lambda, 'p_lambda'=p_lambda))
}

lambda_list <- c(0.2, 1, 1.5, 2)
n_lambda <- length(lambda_list)

N <- 300
AAD_lambda_mat <- matrix(nrow=N, ncol=4)
for(group in 1:N){
  cat('group=', group, '\n')
  q_ex <- data2$qdat[,group]
  r_ex <- data2$rdat[,group]
  GC_AAD_lambda <- matrix(0, nrow=n_lambda, ncol=4)
  colnames(GC_AAD_lambda) <- c('AAD', 'CIV', 'lambda', 'p_lambda')
  for(j in 1:n_lambda){
   GC_AAD_lambda[j,] <- select_lambda_Fun(lambda= lambda_list[j], tau=0.05,method='GC')
  } #[to do]tau=0.01 and change the method
  jlambda <- which.min(GC_AAD_lambda[,'CIV'])
  AAD_lambda_mat[group,] <- GC_AAD_lambda[jlambda,]
}
colnames(AAD_lambda_mat) <- c('AAD', 'CIV', 'lambda', 'p_lambda')
mean(AAD_lambda_mat[,1])#mean(AAD)
sd(AAD_lambda_mat[,1])#sd(AAD)
mean(AAD_lambda_mat[,3])#mean(lambda)
mean(AAD_lambda_mat[,4])#mean(p_lambda)


# real data analysis ------------------------------------------------------
da <- read.table("F:\\最近\\毕业论文\\我的论文\\Code\\CNY_EUR.txt",header = T)
dim(da) #1942
ced <- da[,3]

##对收盘价计算对数收益率
r1<- numeric(length(ced)-1) #1941
for(i in 2:length(ced)){
  r1[i-1] <- (log(ced[i])-log(ced[i-1]))
}

npre <- 200
tau <- 0.05 #[to do] tau=0.01
q <- numeric(length(r1)-npre) #1741
n <- length(q)
r <- numeric(length(q))
for(j in 1:length(q)){
  cat('j = ', j, '\n')
  q[j] <- quantile(r1[(j):(npre+j)], 1-tau) #用来产生和qhat比较的q
  r[j] <- r1[npre+j]
}
c(length(r), length(q))


##UC,CC,Ind
Christoffersen1998 <- function(prob, y, var){
  n <- length(y)
  I <- (y<var)
  n1 <- sum(I) 
  n0 <- n-n1
  n00 <- 0; n01 <- 0; n10 <- 0; n11 <- 0;
  
  for(i in 2:length(y)){
    if(I[i-1]==0 & I[i]==0)n00 <- n00+1
    else if(I[i-1]==0 & I[i]==1)n01 <- n01+1
    else if(I[i-1]==1 & I[i]==0)n10 <- n10+1
    else n11 <- n11+1
  }
  
  pi01 <- ifelse((n00+n01)>0, n01/(n00+n01), 0) 
  pi11 <- ifelse((n10+n11)>0, n11/(n10+n11), 0) 
  pi2 <- (n01+n11)/(n00+n10+n01+n11)
  
  t.uncond <- -2*(n0*log(1-prob) + n1*log(prob) - (n0*log(1-(n1/n)) +n1*log(n1/n)))
  
  t.indep <- -2*log(
    (((1-pi2)^((n00+n10)/2) / (1-pi01)^n00 ) * ((1-pi2)^((n00+n10)/2) / pi01^n01))* 
     ((pi2^((n01+n11)/2) / (1-pi11)^n10) * (pi2^((n01+n11)/2) / pi11^n11 )))

  if(is.na(t.indep)) #If to low values for OS
    t.indep <- -2*((n00+n10)*log(1-pi2)+(n01+n11)*log(pi2) -(n00*log(1-pi01) + 
                             n01*log(pi01) +n10*log(1-pi11) + n11*log(pi11)))

  t.cond <- t.uncond+t.indep
  
  list(t.uncond = t.uncond, t.indep=t.indep, t.cond=t.cond, 
       p.uncond = pchisq(t.uncond, 1, lower.tail=FALSE),
       p.indep = pchisq(t.indep, 1, lower.tail=FALSE),
       p.cond = pchisq(t.cond, 2, lower.tail=FALSE) )
}

library(quantreg)
library(GAS)

## Asymmetric Slope
rt1 <- pmax(0, r)
rt2 <- pmax(0, -r)
dat2 <- data.frame(qt=q[-1], qt1 = q[-n], rt1=rt1[-n], rt2=rt2[-n])

ntr <- 1000
rq1 <- rq(qt~., tau=0.05, data=dat2[1:ntr,])
summary(rq1)

qtesthat <- predict(rq1, dat2[(ntr+1):n, ])
length(qtesthat)#741
qtest <- q[(ntr+1):n]#for what??
Christoffersen1998(tau, r[(ntr+1):(n-1)], -qtesthat[1:740])#740??
BackTest <- BacktestVaR(r[(ntr+1):(n-1)], -qtesthat[1:740], 0.05)
BackTest$DQ

## Natural Cubic Regression Spline
knots <- c(0) #why??
Rmat <- ns(r, knots=knots,intercept = T)
summary(Rmat)
mx <- ncol(Rmat)
theta <- 0.3 #why??
Xfun <- function(j=1){
  Xj <- numeric(n)
  for(t in 2:n){
    for(i in 1:(t-1)){
      Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i, j]
    }
  }
  return(Xj)
}

X <- sapply(1:mx, Xfun)
summary(X)

y <- numeric(n)
for(t in 2:n){
  y[t] <- q[t]-theta^(t-1)*q[1]
}

dat <- data.frame(y=y[-1], x=X[-1,])
names(dat) <- c('y', paste0('x', 1:mx))
rq_NCRS <- rq(y~.+0, tau=0.05, data=dat[1:(ntr),])
summary(rq_NCRS)

ytesthat <- predict(rq_NCRS, dat[(ntr+1):(n), ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
qtest <- q[(ntr+1):n]#for what??
Christoffersen1998(tau, r[(ntr+1):(n-1)], -qtesthat[1:740])
BackTest <- BacktestVaR(r[(ntr+1):(n-1)], -qtesthat[1:740], 0.05)
BackTest$DQ

## Linear Regression Spline
knots <- c(-0.006,0) ##why??
Rmat <- bs(r, knots=knots, degree = 1, intercept = T)
summary(Rmat)
mx <- ncol(Rmat)
theta <- 0.78 ##why??
Xfun <- function(j=1){
  Xj <- numeric(n)
  for(t in 2:n){
    for(i in 1:(t-1)){
      Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i, j]
    }
  }
  return(Xj)
}

X <- sapply(1:mx, Xfun)
summary(X)

y <- numeric(n)
for(t in 2:n){
  y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
}

dat <- data.frame(y=y[-1], x=X[-1,])
names(dat) <- c('y', paste0('x', 1:mx))
rq_LR <- rq(y~.+0, tau=0.05, data=dat[1:ntr,])
summary(rq_LR)

ytesthat <- predict(rq_LR, dat[(ntr+1):n, ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
qtest <- q[(ntr+1):n]#for what??
Christoffersen1998(tau, r[(ntr+1):(n-1)], -qtesthat[1:740])
BackTest <- BacktestVaR(r[(ntr+1):(n-1)], -qtesthat[1:740], 0.05)
BackTest$DQ

## Smoothing Spline
theta <- 0.98 ## why??

y <- numeric(n)
for(t in 2:n){
  y[t] <- q[t]-theta^(t-1)*q[1]
}

ntr <- 1000
dat <- data.frame(y=y[2:ntr], x=r[1:(ntr-1)])
lambda <- 2 ##why??
rq_SS <- rqss(y~ qss(x, constraint = 'N', lambda = lambda) , tau= tau, data=dat)
suma <- summary(rq_SS)
yhat <- predict(rq_SS, dat)
ytesthat <- predict(rq_SS, dat[(ntr+1):n, ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
summary(qtesthat)
qtest <- q[(ntr+1):n]#for what??
Christoffersen1998(tau, r[(ntr+1):(n-1)], -qtesthat[1:740])
BackTest <- BacktestVaR(r[(ntr+1):(n-1)], -qtesthat[1:740], 0.05)
BackTest$DQ

## Indirect Garch
IndirectGARCH<-function(betas, train, CAViaR){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  for(i in 2:length(train)){
    CAViaR[i] <- -(beta1 + beta2 * (CAViaR[i-1]^2) + beta3 * (train[i-1]^2))^(1/2)
  }
  #Objective Function
  res<-sum((tau-(train<CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}
IndirectGARCH(betas=rep(1,3), train=r[1:(ntr-1)], CAViaR=q[2:ntr])
fvx <- optim(rep(1,3), function(betas) IndirectGARCH(betas, train=r[1:(ntr-1)], 
                                                     CAViaR=y[2:ntr]))
betahat <- fvx$par
IndirectGARCHForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- -(beta1+beta2*CAViaR[i-1]^2+beta3*data[i-1]^2)^(1/2)
  }
  return(CAViaR)
}

qtesthat <- IndirectGARCHForecast(betas=betahat, data=r[ntr:(n-1)])
summary(qtesthat)
qtest <- q[(ntr+1):n]#for what??
Christoffersen1998(tau, r[(ntr+1):n], qtesthat)
BackTest <- BacktestVaR(r[(ntr+1):n], qtesthat, 0.05)
BackTest$DQ

## ARMA(1,1)-GARCH(1,1)
library(rugarch)
spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                         garchOrder = c(1, 1), 
                                         submodel = NULL, 
                                         external.regressors = NULL, 
                                         variance.targeting = FALSE), 
                   
                   mean.model     = list(armaOrder = c(1, 1), 
                                         external.regressors = NULL, 
                                         distribution.model = "sstd", 
                                         start.pars = list(), 
                                         fixed.pars = list()))

garch <- ugarchfit(spec = spec, data = q[1:ntr], solver.control = list(trace=0))
pred_garch <- ugarchforecast(garch, n.ahead = length(q)-ntr)
VaR <- quantile(pred_garch, 0.05)
qtest <- q[(ntr+1):n]#for what??
summary(pred_ga)
Christoffersen1998(tau, r[(ntr+1):n], -VaR)

GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                     GASPar = list(location = FALSE, scale = TRUE,
                                   shape = FALSE))
Fit = UniGASFit(GASSpec, r[1:ntr])
Forecast = UniGASFor(Fit, Roll = TRUE, out = r[(ntr+1):n])
VaR <- quantile(Forecast, 0.05)
BackTest <- BacktestVaR(r[(ntr+1):n], VaR, 0.05)
BackTest$DQ
