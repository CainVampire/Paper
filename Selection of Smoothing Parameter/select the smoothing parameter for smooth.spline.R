# https://stackoverflow.com/questions/14929268/how-do-i-select-the-smoothing-parameter-for-smooth-spline
# how do I select the smoothing parameter for smooth.spline()?

# question: I know that the smoothing parameter(lambda) is quite important for fitting a smoothing spline, 
# but I did not see any post here regarding how to select a reasonable lambda (spar=?), I was told that spar 
# normally ranges from 0 to 1. Could anyone share your experience when use smooth.spline()? Thanks.

# agstudy provides a visual way to choose spar. I remember what I learned from linear model class (but not exact) is
# to use cross validation to pick "best" spar. Here's a toy example borrowed from agstudy:

x = seq(1:18)
y = c(1:3,5,4,7:3,2*(2:5),rep(10,4))
splineres <- function(spar){
  res <- rep(0, length(x))
  for (i in 1:length(x)){
    mod <- smooth.spline(x[-i], y[-i], spar = spar)
    res[i] <- predict(mod, x[i])$y - y[i]
  }
  return(sum(res^2))
}

spars <- seq(0, 1.5, by = 0.001)
ss <- rep(0, length(spars))
for (i in 1:length(spars)){
  ss[i] <- splineres(spars[i])
}
plot(spars, ss, 'l', xlab = 'spar', ylab = 'Cross Validation Residual Sum of Squares' , main = 'CV RSS vs Spar')
spars[which.min(ss)]


# Code is not neatest, but easy for you to understand. Also, if you specify cv=T in smooth.spline:
xyspline <- smooth.spline(x, y, cv=T)
xyspline$spar

