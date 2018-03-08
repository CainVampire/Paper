asymetricSlope<-function(betas){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  for(i in 2:length(train)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*max(train[i-1],0)+beta4*max(-train[i-1],0)
  }
  #Objective Function
  res<-sum((tau - (train < CAViaR)) * (train-CAViaR)) / length(train) 
  if(is.na(res)|is.infinite(res)) res<- 1e+10
  #Objective Function
  return(res)
}

asymetricSlopeForecast<-function(betas,data){
  beta1<-betas[1]
  beta2<-betas[2]
  beta3<-betas[3]
  beta4<-betas[4]
  #Create the CAViaR vector
  var<-as.numeric(quantile(data, probs = tau))
  CAViaR<-rep(var,length(data))
  for(i in 2:length(data)){
    CAViaR[i] <- beta1+beta2*CAViaR[i-1]+beta3*max(data[i-1],0)+beta4*max(-data[i-1],0)
  }
  return(CAViaR)
}

