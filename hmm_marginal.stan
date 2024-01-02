functions {

  //Gamma matrix function
  
    matrix G(real G12, real G21) {
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;
    return m;
  }
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

  matrix[nstate, ndept * time] log_omega; // Flatten log_omega for each spatial location and time point

  for (i in 1:ndept) {
    for (t in 1:time) {
      // Assuming spatial dependencies are modeled through lambda
      log_omega[1, (i-1)*time + t] = poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
      log_omega[2, (i-1)*time + t] = poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 1));
    }
  }
}

model {
  // Priors
  G12 ~ beta(1, 1);
  G21 ~ beta(1, 1);
  
  // Likelihood
  
  target += hmm_marginal(log_omega, G(G12, G21), init_density);
}
