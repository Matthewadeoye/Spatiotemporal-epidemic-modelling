library(cmdstanr, quietly = T)
source("Simulation.R")

y<- SimulatedData
init_density<- c(0.4, 0.6)
ndept<- 30
nstate<- 2
e_it<- matrix(rep(1000, ndept*time), ndept, time)
R<- -1*Westmidlands_adjmat
diag(R)<- -rowSums(R, na.rm = T)
qr(R)$rank

#model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/customL.stan"

#model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/hmm_marginal.stan"

model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Spatiotemporal epidemic modelling/customL as function.stan"

model <- cmdstan_model(stan_file = model_code, compile = TRUE)

fit <- model$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, r=r, s=s, u=u, 
                    init_density=init_density, e_it=e_it, R=R), thin = 500, adapt_delta = 0.95, 
                    chains = 4, iter_warmup = 3000, iter_sampling = 2500, parallel_chains = 4) 
fit$summary()
mcmc_trace(fit$draws(variables = c("G12", "G21", "kappa")))
bayesplot::mcmc_scatter(fit$draws(c("G12", "G21", "kappa")), alpha = 0.3)

bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fit$draws(c("G12", "G21", "kappa")))
bayesplot::mcmc_scatter(fit$draws(c("G12", "G21")), alpha = 0.3)
stanfit <- rstan::read_stan_csv(fit$output_files())

#############################################################################################

library(rstan, quietly = T)

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

#Converting data to list
y_list <- lapply(1:nrow(y), function(i) as.matrix(y[i, , drop = FALSE]))

G_matrix<- G(0.2, 0.4)

#Loglikelihood from Stan
Stan_Loglikelihood(y_list, r, s, u, G_matrix, init_density, e_it) 

#Loglikelihood from R
loglikelihood(y,r,s,u,G_matrix,init.density,e.it)