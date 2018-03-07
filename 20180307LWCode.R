# clear the existing environment varable.
rm(list=ls())

# generate data -----------------------------------------------------------

require(fGarch)

tau <- 0.05 ##[to do] tau=0.01
q1 <- qsstd(tau, mean = 0, sd = 2, nu = 6.35, xi = 0.92)
beta0 <- c(-0.30, 0.9, -0.033, -0.066)

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


## test
data1 <- simula(250)   ##[to do] n为500和1000的情况也需要做
q_ex <- data1$qdat[,1] ##[to do] 取第1个样本，下文类似都是第1个样本的结果，需要循环300次
r_ex <- data1$rdat[,1]
q_test <- data1$qdat[,2] 
r_test <- data1$rdat[,2]
summary(data1$rdat)


# parameter estimating ----------------------------------------------------
library(splines)

Rmat <- cbind(1, pmax(0, r_ex), pmax(0,-r_ex))##why not use bs() to generate the design matrix?
summary(r_ex)
# knots <- c(0,1,3)
# Rmat <- bs(r_ex, knots=knots, degree = 1, intercept = T)
# Rmatest <- bs(r_test, knots=knots, degree = 1, intercept = T)
summary(Rmat)
mx <- ncol(Rmat) ##mx equals K (the number of Knots)
theta <- 0.9
n <- 250
Xfun <- function(j=1){
  Xj <- numeric(n)
  for(t in 2:n){ 
    for(i in 1:(t-1)){
      Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i, j]
    }
  }
  return(Xj)
}
X <- sapply(1:mx, Xfun) ## where mx=K, 得出一个矩阵结果;（4.1）中"="右边的第一部分
summary(X)

y <- numeric(n)
for(t in 2:n){
  y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
}


library(quantreg)
dat <- data.frame(y=y[-1], x=X[-1,])
names(dat) <- c('y', 'x1', 'x2', 'x3')
rq1 <- rq(y~.+0, tau=0.05, data=dat)
yhat <- rq1$fitted.values
summary(rq1)

qhat <- yhat + theta^(1:(n-1))*q_ex[1]
AAD_LRS <- mean(abs(qhat - q_ex[2:n]))

betahat_LRS <- c(rq1$coefficients[1],theta,rq1$coefficients[-1]) ##此处theta=beta2,根据
                ##模拟给出的分位数函数，知在CAViaR中theta=beta2,但是在LRS中也相等吗？
names(betahat_LRS) <- paste0('beta',1:4)
betahat_LRS
bias_LRS <- betahat_LRS - beta0

# CAViaR of Asymmetric Slop Model --------------------------------------------
rt2 <- pmax(0, -r_ex)
dat2 <- data.frame(qt=q_ex[-1], qt1 = q_ex[-n], rt1=rt1[-n], rt2=rt2[-n]) ##rt1前面无定义？
head(dat2)
cor(dat2$qt, dat2$rt1)

rq2 <- rq(qt~., tau=0.05, data=dat2)
summary(rq2)

qhat_CAViaR <- rq2$fitted.values
AAD_CAViaR <- mean(abs(qhat_CAViaR - dat2$qt))
betahat_CAViaR <- rq2$coefficients
names(betahat_CAViaR) <- paste0('beta',1:4)
betahat_CAViaR
bias_LRS <- betahat_CAViaR - beta0


# select # of knots ------------------------------------------------------------
n <- 1500
data1 <- simula(n) ##[to do] 下面的simulation只用到了第一种连接函数l(x)
group <- 1  ##[to do] group应该为300，需要循环
q_ex <- data1$qdat[,group] ##前面的模拟中有同名数据集?
r_ex <- data1$rdat[,group] ##前面的模拟中有同名数据集?
summary(data1$rdat)

library(splines)
#Rmat <- cbind(1, pmax(0, r_ex), pmax(0,-r_ex))
summary(r_ex)

knots_list <- list(c(0), ## the knot locations should be find in an automatic way!
                   c(0,3),
                   c(-10,0,3),
                   c(-10,0,3,5),
                   c(-6,-3,0,3,6),
                   c(-6,-3,0,2,3,6))
#attr(bs(r_ex,df=4, degree = 1), 'knots') 
knots <- knots_list[[6]] ## 根据答案指定结点数，感觉泛化能力弱??
Rmat <- bs(r_ex, knots=knots, degree = 1, intercept = T) 
summary(Rmat)
mx <- ncol(Rmat)
theta <- 0.9 
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

library(quantreg)
ntr <- 1000
dat <- data.frame(y=y[2:ntr], x=X[2:ntr,])
names(dat) <- c('y', paste0('x', 1:mx))
rq1 <- rq(y~.+0, tau=0.05, data=dat) ##前面已有rq1?
yhat <- rq1$fitted.values
cbind(yhat, y[2:ntr])[1:20,] #查看前20行
summary(rq1)
X_test_dat <- data.frame(X[(ntr+1):n,])
names(X_test_dat) <- paste0('x', 1:mx)
yhatest <- predict(rq1, X_test_dat)
cbind(ytest, yhatest)[1:20,]

qhat <- yhat + theta^(1:(ntr-1))*q_ex[1]
AAD_LRS <- mean(abs(qhat - q_ex[2:ntr]))


#check function
rhoFun <- function(u, tau){
  tau*u*(u>0) - (1-tau)*u*(u<0)
}

GC <- function(ytest, yhatest, tau){
  mean(rhoFun(ytest-yhatest, tau))
}
GC(ytest, yhatest, 0.05) ##[to do] tau=0.01

SIC <- function(y, yhat, tau, nknots, degree){
  N <- length(y)
  R <- nknots + degree + 1 ## the # of basis function
  log(mean(rhoFun(y-yhat, tau))) + 1/2*log(N)/N*R
}
SIC(y[2:ntr], yhat, 0.05, 1, 1)

AIC <- function(y, yhat, tau, nknots, degree){
  N <- length(y)
  R <- nknots + degree + 1
  log(mean(rhoFun(y-yhat, tau))) + 2*(R+1)/(N-(R+2))
}
AIC(y[2:ntr], yhat, 0.05, 1, 1)

select_knots_Fun <- function(knots, tau=0.05, degree=1,method='GC'){
  nknots = length(knots) ## 此处根据答案给定结点数，并且给定其位置，后用标准比较??
  Rmat <- bs(r_ex, knots=knots, degree = degree, intercept = T)
  #summary(Rmat)
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
  #summary(X)
  y <- numeric(n)
  for(t in 2:n){
    y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
  }
  
  library(quantreg)
  ntr <- 1000
  dat <- data.frame(y=y[2:ntr], x=X[2:ntr,])
  names(dat) <- c('y', paste0('x', 1:mx))
  rq1 <- rq(y~.+0, tau= tau, data=dat) ##前已有同名函数
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
# test
GC1 <- select_knots_Fun(knots=knots_list[[1]], tau=0.05, degree=1,method='GC')
GC4 <- select_knots_Fun(knots=knots_list[[4]], tau=0.05, degree=1,method='GC')

#AAD
nkn_list <- length(knots_list)
GC_AAD <- matrix(0, nrow=nkn_list, ncol=2)
for(j in 1:nkn_list){
  cat('j = ', j, '\n')
  GC_AAD[j,] <- select_knots_Fun(knots=knots_list[[j]], tau=0.05, degree=1,method='GC')
}


# Selection of lambda -----------------------------------------------------
GC1 <- function(ytest, yhatest, tau){
  mean(rhoFun(ytest-yhatest, tau))
}
GC1(ytest, yhatest, 0.05) ##[to do] tau=0.01

SIC1 <- function(y, yhat, tau, R){
  N <- length(y)
  log(mean(rhoFun(y-yhat, tau))) + 1/2*log(N)/N*R
}
SIC1(y[2:ntr], yhat, 0.05,6)

AIC1 <- function(y, yhat, tau, R){
  N <- length(y)
  log(mean(rhoFun(y-yhat, tau))) + 2*(R+1)/(N-(R+2))
}
AIC1(y[2:ntr], yhat, 0.05, 6)

select_lambda_Fun <- function(lambda, tau=0.05,method='GC'){
  #lambda <- 2
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
  rq1 <- myrqss(y~ myqss1(x, constraint = 'N', lambda = lambda, theta=theta), 
                tau= tau, data=dat)
  suma <- summary(rq1)
  yhat <- predict(rq1, dat)
  X_test_dat <- data.frame(x=r_ex[(ntr):(n-1)])
  ytest <- y[(ntr+1):n]
  yhatest <- predict(rq1, X_test_dat)
  qhat <- yhat + theta^(1:(ntr-1))*q_ex[1]
  AAD_LRS <- mean(abs(qhat - q_ex[2:ntr]))
  p_lambda <- suma$edf ## R=p_lambda,means 'effective degree of freedom'
  CIV <- switch(method, 
                GC = GC1(ytest, yhatest, tau),
                SIC = SIC1(y[2:ntr], yhat, tau, p_lambda),
                AIC = AIC1(y[2:ntr], yhat, tau, p_lambda))
  cat('evaluation is finished!------------------------- \n')
  return(c('AAD'=AAD_LRS, 'CIV'=CIV))
}

lambda_list <- c(0.2, 1, 1.5, 2) ##lambda 可以有更automatic way?? Notice table 5.2!
n_lambda <- length(lambda_list)
GC_AAD_labda <- matrix(0, nrow=n_lambda, ncol=2)
for(j in 1:n_lambda){
  cat('j=', j, '\n')
  GC_AAD_labda[j,] <- select_lambda_Fun(lambda= lambda_list[j], tau=0.05,method='GC')
}


# real data analysis ------------------------------------------------------
da <- read.table("F:\\最近\\毕业论文\\我的论文\\Code\\CNY_EUR.txt",header = T)
dim(da) #1942
ced <- da[,3]

##对收盘价计算对数收益率
r1<- numeric(length(ced)-1) #1941
for(i in 2:length(ced)){
  r1[i-1] <- (log(ced[i])-log(ced[i-1]))
}

npre <- 200 ## why 200?这一部分用来做什么?
tau <- 0.05 ##[to do] tau=0.01
q <- numeric(length(r1)-npre) #1741
n <- length(qs) ## 'qs' not found!? 'n' means what?
r <- numeric(length(q))
for(j in 1:length(q)){
  cat('j = ', j, '\n')
  q[j] <- quantile(r1[(j):(npre+j)], 1-tau) ##用来产生和qhat比较的q
  r[j] <- r1[npre+j]
}

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
  
  t.uncond <- -2*( n0*log(1-prob) + n1*log(prob) - ( n0*log(1-(n1/n)) +n1*log(n1/n) ) )
  
  t.indep <- -2*log(
    ( ( (1-pi2)^((n00+n10)/2) / (1-pi01)^n00 ) * ( (1-pi2)^((n00+n10)/2) / pi01^n01 ) ) * 
      ( ( pi2^((n01+n11)/2) / (1-pi11)^n10 ) * ( pi2^((n01+n11)/2) / pi11^n11 )   )
  ) #Factorized for avoiding to low values for OS
  
  if(is.na(t.indep)) #If to low values for OS
    t.indep <- -2*( (n00+n10)*log(1-pi2)+(n01+n11)*log(pi2) -
                      ( n00*log(1-pi01) + n01*log(pi01) +
                          n10*log(1-pi11) + n11*log(pi11) ) )
  #Does not handle n00=0 or n11=0 and is therfore not the default formula for t.indep
  
  t.cond <- t.uncond+t.indep
  
  list(t.uncond = t.uncond, t.indep=t.indep, t.cond=t.cond, 
       p.uncond = pchisq(t.uncond, 1, lower.tail=FALSE),
       p.indep = pchisq(t.indep, 1, lower.tail=FALSE),
       p.cond = pchisq(t.cond, 2, lower.tail=FALSE) )
  
}

## Asymmetric Slope
rt1 <- pmax(0, r)
rt2 <- pmax(0, -r)
dat2 <- data.frame(qt=q[-1], qt1 = q[-n], rt1=rt1[-n], rt2=rt2[-n]) ##前有同名?
##rt1,rt2在这个地方是如何定义的? 是否出现? 前有同名?
head(dat2)
cor(dat2$qt, dat2$rt1)

ntr <- 1000
library(quantreg)
rq1 <- rq(qt~., tau=0.05, data=dat2[1:ntr,]) ##前有同名?
summary(rq1)

qtesthat <- predict(rq1, dat2[(ntr+1):n, ])
qtesthat
qtest <- q[(ntr+1):n]
Christoffersen1998(tau, r[(ntr+1):n], -qtest)


## Natural Cubic Regression Spline
knots <- c(0) ##why??
Rmat <- bs(r, knots=knots, degree = 3, intercept = T) ##应该用ns()??
summary(Rmat)
mx <- ncol(Rmat)
theta <- 0.3 ##why??
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

library(quantreg)
dat <- data.frame(y=y[-1], x=X[-1,])
names(dat) <- c('y', paste0('x', 1:mx))
rq_NCRS <- rq(y~.+0, tau=0.05, data=dat[1:ntr,])
summary(rq_NCRS)

ytesthat <- predict(rq_NCRS, dat[(ntr+1):n, ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
qtest <- q[(ntr+1):n]
Christoffersen1998(tau, r[(ntr+1):n], -qtest)

## Linear Regression Spline
knots <- c(-0.006,0) ##why??
Rmat <- bs(r, knots=knots, degree = 1, intercept = T)
summary(Rmat)
mx <- ncol(Rmat)
theta <- 0.98 ##why??
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

library(quantreg)
dat <- data.frame(y=y[-1], x=X[-1,])
names(dat) <- c('y', paste0('x', 1:mx))
rq_LR <- rq(y~.+0, tau=0.05, data=dat[1:ntr,])
summary(rq_LR)

ytesthat <- predict(rq_LR, dat[(ntr+1):n, ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
qtest <- q[(ntr+1):n]
Christoffersen1998(tau, r[(ntr+1):n], -qtest)

## Smoothing Spline
theta <- 0.98 ##why??
y <- numeric(n)
for(t in 2:n){
  y[t] <- q[t]-theta^(t-1)*q[1]
}

library(quantreg)
ntr <- 1000
dat <- data.frame(y=y[2:ntr], x=r[1:(ntr-1)])
lambda <- 2 ##why??
rq_SS <- rqss(y~ qss(x, constraint = 'N', lambda = lambda) , tau= tau, data=dat)
suma <- summary(rq_SS)
yhat <- predict(rq_SS, dat)
ytesthat <- predict(rq_SS, dat[(ntr+1):n, ])
qtesthat <- ytesthat + theta^(ntr:(n-1))*q[1]
qtest <- q[(ntr+1):n]
Christoffersen1998(tau, r[(ntr+1):n], -qtest)


## Indirect Garch ----------------------------------------------------------

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
fvx <- optim(rep(1,3), function(betas) IndirectGARCH(betas, train=r_ex[1:(ntr-1)], 
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
qtest <- q[(ntr+1):n]
Christoffersen1998(tau, r[(ntr+1):n], -qtest)

