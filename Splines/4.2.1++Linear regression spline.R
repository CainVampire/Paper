library(splines)
#K =14

Rmat2 <- cbind(1, pmax(0, r_ex), pmax(0,-r_ex))
head(Rmat2)

theta <- 0.9
n <- 250

# 此处是把式(4.1)中的那一大坨作为X，其余部分作为y
Xfun <- function(j=1){
  Xj <- numeric(n)
  for(t in 2:n){
    for(i in 1:(t-1)){
      Xj[t] <- Xj[t]+theta^(i-1)*Rmat[t-i,j]
    }
  }
  return(Xj)
}

X <- sapply(1:3, Xfun) #此处是把1，2，3分别代入到了xfun函数中。此处为3的原因是单（内部）节点+2个外部结点=3个结点。
summary(X)


y <- numeric(n)
for(t in 2:n)
  y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]

library(quantreg)
dat <- data.frame(y=y[-1], x=X[-1,]) #此处为什么把y的第一个数去掉，为什么把x的第一个数去掉？
rq1 <- rq(y~.+0, tau=0.05, data=dat)
summary(rq1)


# model 2 -----------------------------------------------------------------
rt1 <- pmax(0, r_ex)
rt2 <- pmax(0, -r_ex)
dat2 <- data.frame(qt=q_ex[-1], qt1 = q_ex[-n], rt1=rt1[-n], rt2=rt2[-n]) #为什么把里面的数去掉？
head(dat2)
cor(dat2$qt, dat2$rt1)
rq2 <- rq(qt~., tau=0.05, data=dat2)
summary(rq2)

#
Rmat <- bs(r_ex, degree = 1, knots = 0, intercept = T)
Rmat <- cbind(1, Rmat[,c(1,3)]) #此处为何意？为什么把1，3列提出和常数列1合并，为什么把2列去掉？why?-20180223
head(Rmat)
summary(Rmat)
