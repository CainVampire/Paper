IndirectGARCH<-function(betas){
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
