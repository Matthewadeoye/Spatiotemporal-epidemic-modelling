functions {
   vector matrix_logsumexp(matrix Mat) {
    int NR = rows(Mat);
    vector[NR] res;
    for (i in 1:NR) {
      res[i] = log_sum_exp(Mat[i,]);      
    }
    return res;
  }  //matrix_logsumexp
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  array[ndept, time] int y;              // data matrix
  vector[time] r;                     // Trend component
  vector[time] s;                     // Seasonal component
  vector[ndept] u;                    // Spatial component  
  simplex[nstate] init_density;        // initial distribution of the Markov chain
  matrix[ndept, time] e_it;            // initial Susceptibles
}


parameters {
  real<lower=0, upper=1> G12;          // transition to hyperendemic
  real<lower=0, upper=1> G21;          // transition to endemic
}

transformed parameters {
    matrix[nstate, nstate] gamma = rep_matrix(0, 2, 2);

    // Build the transition matrix
    gamma[1, 1] = 1-G12;
    gamma[1, 2] = G12;
    gamma[2, 1] = G21;
    gamma[2, 2] = 1-G21;
    

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
}

model {
  // Priors
  G12 ~ beta(1, 1);
  G21 ~ beta(1, 1);
  
  // Likelihood
  
  vector[ndept] row_logsumexp = matrix_logsumexp(logforwardprobs);
  target += sum(row_logsumexp);
}
