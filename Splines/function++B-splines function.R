## 此处默认degree=3，边界值值为0，1，截距为false，即无截距项，
## 其中的边界值为0和1，在我的论文中感觉不是0和1哦
bs <- function(x,degree = 3, interior.knots = NULL, intercept = FALSE, Boundary.knots = range(x)){
  if(missing(x))  stop("You must provide x")
  if(degree < 1)  stop("The spline degree must be at least 1")
  Boundary.knots <- sort(Boundary.knots)
  interior.knots.sorted <- NULL
  if(!is.null(interior.knots)) interior.knots.sorted <- sort(interior.knots)
  ## 此处边界点各重复了（degree+1）次，而论文中重复次数亦为degree+1次！order=degree+1
  knots <- c(rep(Boundary.knots[1],(degree + 1)),interior.knots.sorted,rep(Boundary.knots[2],(degree + 1)))
  cat(knots, '\n')
  ## 此处与原文不一样，原文中K表示内点数目，此处K相当于原文中的K+P+1!!
  K <- length(interior.knots) + degree + 1
  ##经常看到matrix中的第一个数被写为0，是什么意思？2，3个数应依次为行、列
  B.mat <- matrix(0,length(x),K)
  ## ??原论文中的控制点control points 如何在此处表示？原文中的j是从-p开始，此处为1开始！
  ## 但不管在那儿，结点的总个数都等于……
  for(j in 1:K) B.mat[,j]<-basis(x,degree,j,knots)
  if (any(x == Boundary.knots[2]))  B.mat[x == Boundary.knots[2],K] <- 1
  if (intercept == FALSE){
    return(B.mat[,-1])
  }else{
    return(B.mat)
  }
}

#test
#bs(1,degree = 2, interior.knots = c(1,2), intercept = FALSE, Boundary.knots = c(0,3))
#splines::bs(1, degree = 2, knots=c(1,2), intercept = F, Boundary.knots = c(0,3))
