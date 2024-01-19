functions {

   // Intrinsic GMRF density
    real IGMRF1(vector u, real kappa_u, matrix R) {
    int n = rows(R);
    return (((n - 1) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, u));
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
  matrix[ndept, ndept] R;              // Structure matrix (IGMRF1)
}

parameters {
  real<lower=0, upper=1> G12;          // transition to hyperendemic
  real<lower=0, upper=1> G21;          // transition to endemic
  real<lower=0> kappa_u;               // spatial precision parameter
}

transformed parameters {
 vector[ndept] cumHMM;
  array[ndept] matrix[nstate, time] AllOmega;  // AllOmega is an array containing "ndept" matrices having "nstate" rows and "time" columns.
    
  for (i in 1:ndept) {
    for (t in 1:time) {
      AllOmega[i, 1, t] = poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
      AllOmega[i, 2, t] = poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 1));
    }
   cumHMM[i] = hmm_marginal(AllOmega[i], G(G12, G21), init_density);
  }
}

model {
  // Priors
  G12 ~ beta(1, 1);
  G21 ~ beta(1, 1);
  kappa_u ~ gamma(1, 0.01);
  
  // Likelihood
  target += sum(cumHMM) + IGMRF1(u, kappa_u, R);
}
