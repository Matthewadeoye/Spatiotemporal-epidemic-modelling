#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * log(PrecisionR) - PrecisionR/2 * Sumres)
}

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * log(z) - z/2 * Sumres)
}

#Intrinsic GMRF density for spatial components (u_i)
logIGMRF1<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

newlogLikelihood<- function(y, e_it, Model, z_it, z_it2, theta){
  Gamma<- abs(theta[1:2])
  if(Gamma[1]>1)Gamma[1] = 2 - Gamma[1]
  if(Gamma[2]>1)Gamma[2] = 2 - Gamma[2]
  G_mat<- matrix(c(1-Gamma[1], Gamma[1], Gamma[2], 1-Gamma[2]), nrow = length(Gamma), byrow = T)
  r<- theta[5+(1:time)]
  s<- theta[5+((time+1):(2*time))]
  u<- theta[5+(((2*time)+1):((2*time)+ndept))]
  
  if(length(theta) > 5+time+time+ndept){
    ARcoeff<- theta[(5+time+time+ndept+1):length(theta)]
  }else{ARcoeff= NA}
  
  ll<- GeneralLoglikelihood(y, r, s, u, G_mat, e_it, ARcoeff, Model,z_it, z_it2)
}

logPriors<- function(theta, R){
  Gamma<- abs(theta[1:2])
  kappaR<- abs(theta[3])
  kappaS<- abs(theta[4])
  kappaU<- abs(theta[5])
  r<- theta[5+(1:time)]
  s<- theta[5+((time+1):(2*time))]
  u<- theta[5+(((2*time)+1):((2*time)+ndept))]
  
  if(length(theta) > 5+time+time+ndept){
    ARcoeff<- theta[(5+time+time+ndept+1):length(theta)]
  }else{ARcoeff<- c(0, 0)}
  
  lp<- sum(dbeta(Gamma, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
       sum(dgamma(c(kappaR, kappaS), shape = 1, rate = c(0.0001, 0.0001), log = TRUE)) +
       dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
       randomwalk2(r, kappaR) +
       seasonalComp(s, kappaS) +
       logIGMRF1(u, kappaU, R)
      # sum(dexp(ARcoeff, rate = rep(1, length(B)), log = TRUE)) 
}
