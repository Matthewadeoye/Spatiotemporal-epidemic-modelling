#HMC for transition matrix, precision parameters, spatial and temporal components 
set.seed(1223)
#set.seed(1234)
#set.seed(189)
#set.seed(6112)
source("Simulation.R")

y<- SimulationResults[[1]]
plot(1:time, y[10,], type="l")
y
#init_density<- c(0.4, 0.6)
init_density<- c(0.6666667, 0.3333333)
ndept<- length(diag(Westmidlands_adjmat))
#ndept<- length(diag(France_adjmat))
nstate<- 2
e_it<- matrix(rep(1000, ndept*time), ndept, time)
R<- -1*Westmidlands_adjmat
#R<- -1*France_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank

#model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/hmm_marginal.stan"

model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/customL as function.stan"

#model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/GeneralLoglikelihood.stan"

model <- cmdstan_model(stan_file = model_code, compile = TRUE)

initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), r = rep(0, time), s = rep(0, time), kappa_u=20, kappa_r=20, kappa_s=20) 

#initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), r = rep(0, time), s = rep(0, time), kappa_u=20, kappa_r=20, kappa_s=20, B = c(0, 0)) 

fit <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
                    init_density=init_density, e_it=e_it, R=R), init = list(initials, initials, initials, initials), 
                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 


#fit <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 0, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 500, iter_sampling = 1500, parallel_chains = 4, set.seed(1234)) 

#fit1 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 1, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 

#fit2 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 2, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 

#fit3 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 3, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 

#fit4 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 4, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 

#fit5 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 5, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, set.seed(1234)) 

#fit6 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 6, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 500, iter_sampling = 1500, parallel_chains = 4, set.seed(1234)) 

#fit7 <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y,
#                    init_density=init_density, e_it=e_it, R=R, Model = 7, adjacencyMatrix = France_adjmat), init = list(initials, initials, initials, initials), 
#                    chains = 4, iter_warmup = 500, iter_sampling = 1500, parallel_chains = 4, set.seed(1234)) 


fit$summary()
bayesplot::color_scheme_set("purple")
StanUs<- fit$draws(variables = "uconstrained")
mcmc_trace(fit$draws(variables = "uconstrained"))
mcmc_trace(fit$draws(variables = "r"))
mcmc_trace(fit$draws(variables = "s"))

stanfit <- rstan::read_stan_csv(fit$output_files())
fittedG12<- mean(rstan::extract(stanfit)[["G12"]])
fittedG21<- mean(rstan::extract(stanfit)[["G21"]])
fittedKappaU<- mean(rstan::extract(stanfit)[["kappa_u"]])
fittedKappaR<- mean(rstan::extract(stanfit)[["kappa_r"]])
fittedKappaS<- mean(rstan::extract(stanfit)[["kappa_s"]])
fittedUs<- colMeans(rstan::extract(stanfit)[["uconstrained"]])
fittedRs<- colMeans(rstan::extract(stanfit)[["r"]])
fittedSs<- colMeans(rstan::extract(stanfit)[["s"]])
stanfit<- as.matrix(stanfit)
mcmc_recover_hist(stanfit[, 100:159], r)
mcmc_recover_hist(stanfit[, 160:219], s)
mcmc_recover_hist(stanfit[, 223:317], u)
mcmc_recover_hist(stanfit[, 1:2], c(0.2,0.4))

mcmc_trace(fit$draws(variables = c("G12", "G21", "kappa_u", "kappa_r", "kappa_s")))
bayesplot::mcmc_dens(fit$draws(c("G12", "G21", "kappa_u")))
bayesplot::mcmc_scatter(fit$draws(c("G12", "G21")), alpha = 0.3)
stanfit <- rstan::read_stan_csv(fit$output_files())

#Data from paper
library(surveillance)
MeningococcalData<- data("meningo.age")
plot(meningo.age, title="Meningococcal infections in France 1985-97")
MeningococcalData<- disProg2sts(meningo.age, map=NULL)
MeningococcalIncidence<- observed(MeningococcalData)
plot(1:156, rowSums(MeningococcalIncidence), type="l")

#Results and figures for paper for general model
stanfit <- rstan::read_stan_csv(fit$output_files())

plot(1:time, colSums(y), type="l", xlab="Months", ylab="Total observed counts")

plot(1:time, fittedRs, type="l", xlab="Months", ylab="Fitted trend components")

plot(1:time, fittedSs, type="l", xlab="Months", ylab="Fitted Seasonal component")

expUi<- exp(fittedUs)

rrDF<- data.frame(rr=expUi, nuts3=France.shp$nuts3, id=c(1:ndept), AreaName=France.shp$nom, months=c(1:ndept))

shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/ModifiedFranceShapefile.shp")
shapefile <- fortify(shapefile, region = 'nuts3')

shp_rrDF<- sp::merge(shapefile, rrDF,
                               by.x="id",
                               by.y="nuts3",
                               all.x=F)
shp_rrDF <- arrange(shp_rrDF, order)

rr.map <- ggplot(shp_rrDF, aes(long, lat, fill = rr, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.3, end=1) +  # Reverse the color scale
  coord_equal() + theme_void() +
  # ggtitle('') +
  labs(fill = "Median relative risk")
rr.map


source("OutbreakProbability.R")
image(SimulationResults[[2]])
image(Decoding(y, fittedRs, fittedSs, fittedUs, G(fittedG12, fittedG21), init_density, e.it))

#ROC curve
library(pROC)
true.outbreaks <- as.vector(SimulationResults[[2]])
predicted.outbreaks <- as.vector(Decoding(y, fittedRs, fittedSs, fittedUs, G(fittedG12, fittedG21), init_density, e.it))
roc_data <- roc(true.outbreaks, predicted.outbreaks)
plot(roc_data, main = "ROC Curve", col = "blue", lwd = 2)
legend("bottomright", legend = paste("AUC =", round(auc(roc_data), 2)), col = "blue", lty = 1)


#############################################################################################

expose_stan_functions(stanc(model_code = "functions {

  //Gamma matrix function
  
    matrix G(real G12, real G21) {
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;
    return m;
  }


  //loglikelihood via forward filtering

  real Stan_Loglikelihood(array[,] int y, vector r, vector s, vector u, matrix gamma, vector init_density, matrix e_it) {

  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

  matrix[ndept, nstate] logforwardprobs;
  matrix[time, nstate] Alphas;
  matrix[time, nstate] Betas;
  matrix[time, nstate] alpha;
  matrix[time, nstate] beta;

  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 1));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 1));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 1));
         Alphas[t, 1] = log_sum_exp(alpha[t,]);
         Alphas[t, 2] = log_sum_exp(beta[t,]);
}
      
      // save the final forward probabilities

      logforwardprobs[i,1] = log_sum_exp(alpha[time,]);
      logforwardprobs[i,2] = log_sum_exp(beta[time,]);

    }

      vector[ndept] rowlogsumexp;
      for (i in 1:ndept) {
      rowlogsumexp[i] = log_sum_exp(logforwardprobs[i,]);      
    }
      real fullLogLikelihood = sum(rowlogsumexp);
      return fullLogLikelihood;
  }
}"))

#G function from Stan
G(0.2, 0.4)
G_matrix<- G(0.2, 0.4)

#Converting matrix y to lists of matrices
y_list <- lapply(1:nrow(y), function(i) as.matrix(y[i, , drop = FALSE]))

adj_list <- lapply(1:nrow(France_adjmat), function(i) as.matrix(France_adjmat[i, , drop = FALSE]))

#Loglikelihood from Stan
Stan_Loglikelihood(y_list, r, s, u, G_matrix, init_density, e_it) 

Stan_Loglikelihood(y_list, r, s, u, G_matrix, init_density, e_it, c(0.1, 0.2), 0, adj_list) 

#Loglikelihood from R
source("GeneralLoglikelihood.R")
loglikelihood(y,r,s,u,G_matrix,init.density,e.it)
GeneralLoglikelihood(y,r,s,u,G_matrix,init.density,e.it,B=c(0.1, 0.2),model = 0,adjacencyMatrix= France_adjmat)

expose_stan_functions(stanc(model_code = "functions {
 real hmarginal (matrix logOmega, matrix GammaMatrix, vector InitialDensity){
 return hmm_marginal(logOmega, GammaMatrix, InitialDensity);
 }
}"))

cumHMM<- 0
for (i in 1:ndept) {
  AllOmega = matrix(data= NA, nrow=2, ncol=time)
  for (t in 1:time) {
    AllOmega[1, t] = dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 0), log=T)
    AllOmega[2, t] = dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + 1), log=T)
  }
  cumHMM = cumHMM + hmarginal(AllOmega, G_matrix, init_density)
}
cumHMM