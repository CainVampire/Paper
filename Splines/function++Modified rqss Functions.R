#This code modifies the quantile smoothing spline functions in quantreg package 
#to estimate the recursive quantile regression model when u=1 and v=1.

library(quantreg)
library(SparseM)

# Given a time series data {r_t}, y is the current obs and x is the lag-1 day obs

make_new = function(x,y){
  #make "y" a new "y" to be used as the response
  sx = sort(x)
  loc = which(sx[-length(sx)]==sx[-1])
  add_len = length(loc)
  c(y,rep(0,add_len))
}


# myqss1 generates the design matrix for rqss fit, to be used in myrqss

myqss1 <- function(x,constraint="N",lambda=1,dummies=dummies,theta=theta,ndum=0,w=rep(1,length(x))){
  xun <- unique(x[order(x)])
  h <- diff(xun)
  nh <- length(h)
  nx <- length(x)
  p <- nh + 1
  B <- new("matrix.csr",ra = c(rbind(-1/h,1/h)))
  ja = as.integer(c(rbind(1:nh, 2:(nh + 1))))
  ia = as.integer(2 * (1:(nh + 1)) - 1)
  dimension = as.integer(c(nh, nh + 1))
  makeD <- function(p) {
    new("matrix.csr", ra = c(rbind(rep(-1, (p - 1)), rep(1,(p - 1)))),
    ja = as.integer(c(rbind(1:(p - 1), 2:p))),
    ia = as.integer(2 * (1:p) - 1),
    dimension = as.integer(c(p - 1, p)))
  }
  D <- makeD(nh)  
  A <- D %*% B
  # up to this point, the design matrix is the same as in the original rqss
  
  if (length(xun)<nx){
    # this part adjusts the design matrix for identical values in x
    num_same = function(u){sum(sx==u)}
    num = sapply(xun,num_same)
    loc_same = c(which(num >= 2),0)
    Anew = as.matrix(A[,1])
    for (k in 1:ncol(A)) {
      Anew = cbind(Anew,matrix(rep(as.matrix(A[,k]/num[k]),num[k]), ncol=num[k]))
    }
    A=Anew[,-1]
  }
  ox = order(x)   #this part does the column manipulation in the
  rankx = 1       #design matrix for estimating recursive regression
  rankx[ox] = 1:nx
  A <- A[,rankx]
  A <- cbind( A[,-nx]-theta*A[,-1] , A[,nx])
  A <- as.matrix.csr(A)
  design = diag(1,nrow=nx)
  if (length(xun)<nx) {
    loc = which(sx[-length(sx)]==sx[-1])
    add = matrix(0, nrow=length(loc), ncol=nx)
    for (i in 1:length(loc)){
      add[i,loc[i]] = 1
      add[i,loc[i]+1] = -1
    }
    add <- add[,rankx]
    if (length(loc)==1) {
      add <- c( add[-nx]-theta*add[-1] , add[nx])
    } else {
      add <- cbind( add[,-nx]-theta*add[,-1] , add[,nx])
    }
    design = rbind(as.matrix(design), add)
  }
  
  F<-as.matrix.csr(design)
  
  switch (constraint,V = {
    R <- A
    r <- rep(0,nrow(R))
  }, C = {
    R <- -A
    r <- rep(0,nrow(R))
  }, I = {
    R <- makeD(p)
    r <- rep(0, p - 1)
  }, D = {
    R <- -makeD(p)
    r <- rep(0, p - 1)
  }, VI = {
    R <- makeD(p)
    R <- rbind(R,A)
    r <- rep(0, nrow(R))
  },VD = {
    R <- -makeD(p)
    R <- rbind(R,A)
    r <- rep(0, nrow(R))  
  }, CI = {
    R <- makeD(p)
    R <- rbind(R,-A)
    r <- rep(0, nrow(R))
  }, CD = {
    R <- -makeD(p)
    R <- rbind(R,-A)
    r <- rep(0, nrow(R))  
  }, N = {
    R = NULL
    r = NULL
    }
  )
  list(x = list(x = xun), F = F, lambda = lambda, A = A, R = R[,-1], r = r)
}


qss<-function (x, constraint = "N", lambda = 1, ndum = 0, dummies = NULL,theta,w = rep(1, length(x))){
  if (is.matrix(x)) {
    if (ncol(x) == 2)
      qss <- qss2(x, constraint = constraint, dummies = dummies, lambda = lambda, ndum = ndum, w = w)
    else if (ncol(x) == 1)
      x <- as.vector(x)
    else stop("qss objects must have dimension 1 or 2")
  }
  if (is.vector(x))
    qss <- myqss1(x, constraint = constraint, lambda = lambda, theta=theta, dummies = dummies, ndum = ndum, w = w)
  qss
}


#only changed the part of rqss fitting on qss terms without restrictions

myrqss<-function (formula, tau = 0.5, data = parent.frame(), weights, na.action, method = "fn", contrasts = NULL, ...){
  call <- match.call()
  m <- match.call(expand = FALSE)
  temp <- c("", "formula", "data", "weights", "na.action")
  m <- m[match(temp, names(m), nomatch = 0)]
  m[[1]] <- as.name("model.frame")
  special <- "qss"
  Terms <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  qssterms <- attr(Terms, "specials")$qss
  dropx <- NULL
  
  
  if (length(qssterms)) {
    tmpc <- untangle.specials(Terms, "qss")
    ord <- attr(Terms, "order")[tmpc$terms]
    if (any(ord > 1))
      stop("qss can not be used in an interaction")
    dropx <- tmpc$terms
    if (length(dropx))
      Terms <- Terms[-dropx]
    attr(Terms, "specials") <- tmpc$vars
    qssnames <- unlist(lapply(parse(text = tmpc$vars), function(x) deparse(x[[2]])))
  }
  
  m$formula <- Terms
  m <- eval(m, parent.frame())
  weights <- model.extract(m, weights)
  process <- (tau < 0 || tau > 1)
  Y <- model.extract(m, "response")
  X <- model.matrix(Terms, m, contrasts)
  p <- ncol(X)
  pf <- environment(formula)
  
  if (length(qssterms) > 0){
    F <- as.matrix.csr(X)
    qss <- lapply(tmpc$vars, function(u) eval(parse(text = u), data, enclos = pf))
    mqss <- length(qss)
    ncA <- rep(0, mqss + 1)
    nrA <- rep(0, mqss + 1)
    nrR <- rep(0, mqss + 1)
    for (i in 1:mqss) {
      F <- cbind(F, qss[[i]]$F)
      ncA[i + 1] <- ncol(qss[[i]]$A)
      nrA[i + 1] <- nrow(qss[[i]]$A)
      nrR[i + 1] <- ifelse(is.null(nrow(qss[[i]]$R)), 0, nrow(qss[[i]]$R))
    }
    
    F = F[, -1]
    A <- as.matrix.csr(0, sum(nrA), sum(ncA))
    if (sum(nrR) > 0) {
      R <- as.matrix.csr(0, sum(nrR), sum(ncA))
      nrR <- cumsum(nrR)
    }
    
    ncA <- cumsum(ncA)
    nrA <- cumsum(nrA)
    lambdas <- rep(0, mqss)
    for (i in 1:mqss) {
      lambdas[i] <- qss[[i]]$lambda
      Arows <- (1 + nrA[i]):nrA[i + 1]
      Acols <- (1 + ncA[i]):ncA[i + 1]
      A[Arows, Acols] <- qss[[i]]$lambda * qss[[i]]$A
      if (nrR[i] < nrR[i + 1])
        R[(1 + nrR[i]):nrR[i + 1], (1 + ncA[i]):ncA[i + 1]]<-qss[[i]]$R
    }
    
    if (nrR[mqss + 1] > 0) {
      R <- cbind(as.matrix.csr(0, nrR[mqss + 1], p), R)
      r <- rep(0, nrR[mqss + 1])
    }
    else {
      R <- NULL
      r <- NULL
    }
    
    X <- rbind(F, A)
    Y <- c(Y, rep(0, nrow(A)))
    rhs <- t(rbind((1 - tau) * F, 0.5 * A)) %*% rep(1, nrow(X))
    XpX <- t(X) %*% X
    nnzdmax <- XpX@ia[length(XpX@ia)] - 1
    nsubmax <- max(nnzdmax, floor(1000 + exp(-1.6) * nnzdmax^1.2))
    nnzlmax <- floor(2e+05 - 2.8 * nnzdmax + 7e-04 * nnzdmax^2)
    tmpmax <- floor(1e+05 + exp(-12.1) * nnzdmax^2.35)
    fit <- if (length(r) > 0)
      rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfnc", R = R, r = r, nsubmax = nsubmax, nnzlmax = nnzlmax, tmpmax = tmpmax)
    else rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfn", nnzlmax= nnzlmax,nsubmax = nsubmax, tmpmax = tmpmax)
    
    for (i in 1:mqss) {
      ML <- p + ncA[i]
      MU <- p + ncA[i + 1] - 1
      qss[[i]] <- list(xyz = cbind(qss[[i]]$x$x, qss[[i]]$x$y), coef = fit$coef[ML:MU], dummies = qss[[i]]$dummies, X=as.matrix(F))
    }
    names(qss) <- qssnames
    fit$qss <- qss
  }
  else{
    fit <- if (length(weights))
      rq.wfit(X, Y, tau = tau, weights, method, ...)
    else rq.fit(X, Y, tau = tau, method, ...)
  }
  
  fit$terms <- Terms
  fit$formula <- formula
  fit$tau <- tau
  if (length(qssterms)) {
    fit$lambdas <- lambdas
    fit$nrA <- nrA
  }
  else fit$lambdas <- fit$nrA <- NA
  attr(fit, "na.message") <- attr(m, "na.message")
  fit
}

