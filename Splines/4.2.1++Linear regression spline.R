library(splines)
#K =14
Rmat2 <- cbind(1, pmax(0, r_ex), pmax(0,-r_ex))
head(Rmat2)
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
X <- sapply(1:3, Xfun)
summary(X)
y <- numeric(n)
for(t in 2:n)
  y[t] <- q_ex[t]-theta^(t-1)*q_ex[1]
library(quantreg)
dat <- data.frame(y=y[-1], x=X[-1,])
rq1 <- rq(y~.+0, tau=0.05, data=dat)
summary(rq1)


# model 2 -----------------------------------------------------------------
rt1 <- pmax(0, r_ex)
rt2 <- pmax(0, -r_ex)
dat2 <- data.frame(qt=q_ex[-1], qt1 = q_ex[-n], rt1=rt1[-n], rt2=rt2[-n])
head(dat2)
cor(dat2$qt, dat2$rt1)
rq2 <- rq(qt~., tau=0.05, data=dat2)
summary(rq2)

#
Rmat <- bs(r_ex, degree = 1, knots = 0, intercept = T)
Rmat <- cbind(1, Rmat[,c(1,3)]) #此处为何意？why?-20180223
head(Rmat)
summary(Rmat)
