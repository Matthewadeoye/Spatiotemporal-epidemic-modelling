source("GeneralLoglikelihood.R")
source("logDensities.R")

ModelEvidence<- function(y, e_it, R, Model, z_it, z_it2, posteriorSamples, num_samples = 10000){
  
  if(Model == 0){
    posteriorSamples<- posteriorSamples[,-((ncol(posteriorSamples)-1): ncol(posteriorSamples))]
  }else if(Model == 1 || Model == 2 || Model == 4 || Model ==5){
    posteriorSamples<- posteriorSamples[,-ncol(posteriorSamples)] 
  }else if(Model == 3 || Model ==6){
    posteriorSamples<- posteriorSamples
  }
  
  mu<- colMeans(posteriorSamples)
  varcov<- cov(posteriorSamples)
  
  EvaluatedLogDensities<- numeric(num_samples)

for(i in 1:num_samples){
  scaled_sigma<- varcov * (length(mu) - 2)/length(mu)
  theta<- mvtnorm::rmvt(n=1, sigma = scaled_sigma, df = length(mu), delta = mu)

#Evaluate LogDensities
  EvaluatedLogDensities[i]<- newlogLikelihood(y, e_it, Model, z_it, z_it2, theta) + logPriors(theta, R)
                            - sum(mvtnorm::dmvt(theta, delta = mu, sigma = scaled_sigma, df = length(mu), log = TRUE)) 
}
  MarginalLikelihood<- log(1) - log(num_samples) + logSumExp(EvaluatedLogDensities)
  return(MarginalLikelihood)
}