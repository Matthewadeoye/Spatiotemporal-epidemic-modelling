functions {

  // Intrinsic GMRF density
    real IGMRF1(vector u, real kappa_u, matrix R) {
    int n = rows(R);
    return (((n - 1) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, u));
  }

  // or

    real IGMRF2(vector u, matrix Q) {
    int n = rows(Q);
    vector[n] EVs = eigenvalues_sym(Q);
    real product_nonzero_EVs = 1.0;
    for (i in 1:n) {
    if (EVs[i] >= 0.00000001) {
      product_nonzero_EVs *= EVs[i];
    }
  }
    return (-(((n - 1) / 2.0) * log(2.0 * pi())) + (0.5 * log(product_nonzero_EVs)) - (0.5 * quad_form(Q, u)));
  }

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
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  array[ndept, time] int y;           // data matrix
  vector[time] r;                     // Trend component
  vector[time] s;                     // Seasonal component
  vector[ndept] u;                    // Spatial component  
  simplex[nstate] init_density;       // initial distribution of the Markov chain
  matrix[ndept, time] e_it;           // initial Susceptibles
  matrix[ndept, ndept] R;             // Structure / Precision matrix (IGMRF1/IGMRF2)
}

parameters {
  real<lower=0, upper=1> G12;          // transition to hyperendemic
  real<lower=0, upper=1> G21;          // transition to endemic
  real<lower=0> kappa_u;                 // precision parameter
  // vector[ndept] u;                    // Spatial component 
}

model {
  // Priors
  G12 ~ beta(1, 1);
  G21 ~ beta(1, 1);
  kappa_u ~ gamma(1, 0.01);     

  // Likelihood
  target += Stan_Loglikelihood(y, r, s, u, G(G12, G21), init_density, e_it) + IGMRF1(u, kappa_u, R);
}
