#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#check_cmdstan_toolchain()
#install_cmdstan(cores = 2)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
cmdstan_path()

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
mod$print()

# names correspond to the data block in the Stan program
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

fit$summary()
fit$summary(variables = c("theta", "lp__"), "mean", "sd")

# use a formula to summarize arbitrary functions, e.g. Pr(theta <= 0.5)
fit$summary("theta", pr_lt_half = ~ mean(. <= 0.5))

# summarise all variables with default and additional summary measures
fit$summary(
  variables = NULL,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)

# default is a 3-D draws_array object from the posterior package
# iterations x chains x variables
draws_arr <- fit$draws() # or format="array"
str(draws_arr)

# draws x variables data frame
draws_df <- fit$draws(format = "df")
str(draws_df)
print(draws_df)

# this should be identical to draws_df created via draws(format = "df")
draws_df_2 <- as_draws_df(draws_arr)
identical(draws_df, draws_df_2)

mcmc_hist(fit$draws("theta"))


fit_mle <- mod$optimize(data = data_list, seed = 123)
fit_mle$print() # includes lp__ (log prob calculated by Stan program)

mcmc_hist(fit$draws("theta")) +
  vline_at(fit_mle$mle("theta"), size = 1.5)


############################################################################################
#SIR model - Rstan package

library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
# Automatically save compiled Stan models so they can be ran multiple
# times without getting recompiled :

rstan_options(auto_write = TRUE)

# Chains will run in parallel when possible :
options(mc.cores = parallel :: detectCores())

onset<- influenza_england_1978_school$date # Onset date

cases<- influenza_england_1978_school$in_bed # Number of infected students

N=length(onset) # Number of days observed throughout the outbreak

pop=763 # Population

sample_time=1: N

# Modify data into a form suitable for Stan
flu_data = list( n_obs = N ,
                 n_theta = 2,
                 n_difeq = 3,
                 n_pop = pop ,
                 y = cases ,
                 t0 = 0,
                 ts = sample_time )

# Specify parameters to monitor
parameters =c("y_hat", "y_init", "theta" , "R_0")


n_chains =4
n_warmups =500
n_iter =100500
n_thin =50
set.seed(1234)

# Set initial values :
ini = function(){
  list(theta =c(runif(1 ,0 ,5),runif(1,0.2,0.4)) ,
       S0 = runif(1,(pop-3)/pop,(pop-1)/pop))
}

nuts_fit = stan(file = "C:/Users/Matthew Adeoye/Documents/SIR.stan" , # Stan program
                data = flu_data , # list of data
                pars = parameters , # monitored parameters
                init = ini , # initial parameter values
                chains = n_chains , # number of Markov chains to run
                warmup = n_warmups , # number of warmup iterations per chain
                iter = n_iter , # number of iterations per chain (+ warmup )
                thin = n_thin , # period for saving samples
                seed =13219)


print(nuts_fit)
nuts_fit_summary <- summary(nuts_fit, pars = c("lp__", "theta[1]","theta[2]","y_init[1]", "R_0"))$summary
print(nuts_fit_summary, scientific = FALSE, digits =2)
# Obtain the generated samples :
posts<- rstan::extract(nuts_fit)

library(bayesplot)
posterior<- as.array(nuts_fit)
mcmc_trace(posterior, pars =c("lp__", "theta[1]", "theta[2]", "y_init[1]", "R_0"))
pairs(nuts_fit, pars= c("theta[1]", "theta[2]", "y_init[1]"), labels= c("beta", "gamma", "s(0)"),cex.labels=1.5,font.labels=9,
      condition= "accept_stat__")
