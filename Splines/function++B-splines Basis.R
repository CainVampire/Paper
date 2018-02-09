## recursion formula is used to construct B-splines.
basis <- function(x,degree,i,knots){
  if (degree == 0){
    B <- ifelse((x >= knots[i]) & (x<knots[i+1]),1,0)
  }else{
    if (knots[degree+i]-knots[i] == 0){
      temp1 <- 0
    }else{
      temp1 <- (x-knots[i])/(knots[degree+i]-knots[i])
    }
    if (knots[i+degree+1]-knots[i+1] == 0){
      temp2 <- 0
    }else{
      temp2 <- (knots[i+degree+1]-x)/(knots[i+degree+1]-knots[i+1])
    }
    B <- temp1*basis(x,(degree-1),i,knots)+temp2*basis(x,(degree-1),(i+1),knots)
  }
  return(B)
}
