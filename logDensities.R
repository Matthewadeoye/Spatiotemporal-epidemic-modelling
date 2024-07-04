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
  if(Model==0){
    Gamma<- c(0.5, 0.5)
    G_mat<- matrix(c(1-Gamma[1], Gamma[1], Gamma[2], 1-Gamma[2]), nrow = length(Gamma), byrow = T)
    r<- theta[3+(1:time)]
    s<- theta[3+((time+1):(2*time))]
    u<- theta[3+(((2*time)+1):((2*time)+ndept))]
    ARcoeff<- c(0, 0)
  }else{
    Gamma<- theta[1:2] 
    G_mat<- matrix(c(1-Gamma[1], Gamma[1], Gamma[2], 1-Gamma[2]), nrow = length(Gamma), byrow = T)
    r<- theta[5+(1:time)]
    s<- theta[5+((time+1):(2*time))]
    u<- theta[5+(((2*time)+1):((2*time)+ndept))]
    ARcoeff<- theta[(5+time+time+ndept+1):length(theta)]
  }
  
  ll<- GeneralLoglikelihood_cpp(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
  return(ll)
}

logPriors<- function(theta, R, Model){
  if(Model==0){
    kappaR<- theta[1]
    kappaS<- theta[2]
    kappaU<- theta[3]
    r<- theta[3+(1:time)]
    s<- theta[3+((time+1):(2*time))]
    u<- theta[3+(((2*time)+1):((2*time)+ndept))]
    lp<- sum(dgamma(c(kappaR, kappaS), shape = 1, rate = c(0.0001, 0.0001), log = TRUE)) +
         dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
         randomwalk2(r, kappaR) +
         seasonalComp(s, kappaS) +
         logIGMRF1(u, kappaU, R)
  }else{
    Gamma<- theta[1:2]
    kappaR<- theta[3]
    kappaS<- theta[4]
    kappaU<- theta[5]
    r<- theta[5+(1:time)]
    s<- theta[5+((time+1):(2*time))]
    u<- theta[5+(((2*time)+1):((2*time)+ndept))]
    ARcoeff<- theta[(5+time+time+ndept+1):length(theta)] 
    lp<- sum(dbeta(Gamma, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
         sum(dgamma(c(kappaR, kappaS), shape = 1, rate = c(0.0001, 0.0001), log = TRUE)) +
         dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
         randomwalk2(r, kappaR) +
         seasonalComp(s, kappaS) +
         logIGMRF1(u, kappaU, R) + 
         sum(dexp(ARcoeff, rate = rep(1, length(ARcoeff)), log = TRUE))
    }
  return(lp)
}
