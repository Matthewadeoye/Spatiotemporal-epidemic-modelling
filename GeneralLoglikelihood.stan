functions {
  // Neighbors
  array[] int neighbors(int current, array[,] int adjacencyMatrix) {
  int ndept = dims(adjacencyMatrix)[1];
  array[ndept] int indexes;
  int num_neighbors = 0;
  
      for (a in 1:ndept) {
        if (current != a && adjacencyMatrix[current, a] > 0) {
          indexes[a] = a;
          num_neighbors += 1;
        } else {
          indexes[a] = 0;
        }
      }
  array[num_neighbors] int cleanedIndexes;
  int g = 1;
  for (j in 1:num_neighbors){
  for (i in g:ndept) {
    g += 1;
    if (indexes[i] > 0) {
      cleanedIndexes[j] = indexes[i];
      break;
     }
   }
 }
    return cleanedIndexes;
}
  
  // Intrinsic GMRF density
  real IGMRF1_lpdf(vector uconstrained, real kappa_u, matrix R) {    
    int n = rows(R);
    return (((n - 1) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, uconstrained));
  }
  
  // Random walk density
  real randomwalk2_lpdf(vector r, real kappa_r) { 
    int time = dims(r)[1];  
    real res = 0;
    for (i in 3:time) {
      res += (r[i-2] - 2 * r[i-1] + r[i])^2;
    }
    return (((time - 2) / 2.0) * log(kappa_r) - (kappa_r / 2.0) * res);   
  }
  
  // Seasonal components' density
    real seasonalComp_lpdf(vector s, real kappa_s) { 
    int time = dims(s)[1];  
    real res = 0;
    for (i in 12:time) {
    res += (sum(s[(i-11):(i-0)]))^2;
  }
    return (((time - 11) / 2.0) * log(kappa_s) - (kappa_s / 2.0) * res);   
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

  real Stan_Loglikelihood(array[,] int y, vector r, vector s, vector u, matrix gamma, vector init_density, matrix e_it, vector B, int Model, array[,] int adjacencyMatrix) {

  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstate = rows(gamma);

  matrix[ndept, nstate] logforwardprobs;
  matrix[time, nstate] Alphas;
  matrix[time, nstate] Betas;
  matrix[time, nstate] alpha;
  matrix[time, nstate] beta;
  
  if(Model == 0){

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
  
  else if(Model == 1){
  //Model1
    matrix[ndept, time] z_it;

    for(i in 1:ndept){
        z_it[i, 1] = 0;
      for(t in 2:time){
        if(y[i, t-1] > 0){
          z_it[i, t] = 1;
        } else {
          z_it[i, t] = 0;
        }
      }
    }
  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1]));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t]));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t]));
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
  
  else if(Model == 2){
  //Model2
     matrix[ndept, time] z_it;
     array[ndept] int length;   
      
      for(i in 1:ndept){
      length[i] = dims(neighbors(i, adjacencyMatrix))[1];
            z_it[i, 1] = 0;
      
        for(t in 2:time){
          for(b in 1:length[i]){
          if(y[i, t-1] > 0 || y[neighbors(i, adjacencyMatrix)[b], t-1] > 0){
            z_it[i, t]= 1;
          } else {
            z_it[i, t]= 0;
          }
        }
      }
    }
      
  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1]));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t]));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t]));
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
  
  else if(Model == 3){
  //Model3
     matrix[ndept, time] z_it;
     for(i in 1:ndept){
          z_it[i, 1] = 0;
      for(t in 2:time){
        if(y[i, t-1] > 0){
          z_it[i, t] = 1;
        } else {
          z_it[i, t] = 0;
        }
      }
    }
     matrix[ndept, time] z_it2;
     array[ndept] int length;   
      
      for(i in 1:ndept){
      length[i] = dims(neighbors(i, adjacencyMatrix))[1];
           z_it2[i, 1] = 0;
      
        for(t in 2:time){
          for(b in 1:length[i]){
          if(y[neighbors(i, adjacencyMatrix)[b], t-1] > 0){
            z_it2[i, t]= 1;
          } else {
            z_it2[i, t]= 0;
          }
        }
      }
    }
      
  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1] + B[2] * z_it2[i, 1]));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t] + B[2] * z_it2[i, t]));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i, t] + B[2] * z_it2[i, t]));
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
    
    else if(Model == 4){

  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 1)));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1)));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1)));
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
  
   else if(Model == 5){
  //Model5
     matrix[ndept, time] z_it;
     array[ndept] int length;   
     array[ndept] int sumNeighbor;
      
      for(i in 1:ndept){
      length[i] = dims(neighbors(i, adjacencyMatrix))[1];
      z_it[i, 1] = 0;
      
        for(t in 2:time){
            sumNeighbor[i] = 0;
          for(b in 1:length[i]){
            sumNeighbor[i] += y[neighbors(i, adjacencyMatrix)[b], t-1];
          }
          z_it[i, t] = sumNeighbor[i];
        }
      }
      
  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 0 + 1)));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + z_it[i, t] + 1)));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + z_it[i, t] + 1)));
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
  
  else if(Model == 6){
  //Model6
     matrix[ndept, time] z_it;
     array[ndept] int length;   
     array[ndept] int sumNeighbor;
      
      for(i in 1:ndept){
      length[i] = dims(neighbors(i, adjacencyMatrix))[1];
      z_it[i, 1] = 0;
      
        for(t in 2:time){
            sumNeighbor[i] = 0;
          for(b in 1:length[i]){
            sumNeighbor[i] += y[neighbors(i, adjacencyMatrix)[b], t-1];
          }
          z_it[i, t] = sumNeighbor[i];
        }
      }
      
  for (i in 1:ndept) {
  
  // Initialization of the first time step for each department

    alpha[1, 1] = log(init_density[1]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + 0));
    alpha[1, 2] = 0;
    beta[1, 1] = 0;
    beta[1, 2] = log(init_density[2]) + poisson_lpmf(y[i, 1] | e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * log(0 + 1) + B[2] * (z_it[i, 1] + 1)));
    
    Alphas[1, 1] = alpha[1, 1];
    Alphas[1, 2] = beta[1, 2];    
  
  // Dynamic programming loop for the remaining time steps
      for (t in 2:time) {

         alpha[t, 1] = Alphas[t-1, 1] + log(gamma[1, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         alpha[t, 2] = Alphas[t-1, 2] + log(gamma[2, 1]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + 0));
         beta[t, 1] =  Alphas[t-1, 1] + log(gamma[1, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1) + B[2] * (z_it[i, t] + 1)));
         beta[t, 2] =  Alphas[t-1, 2] + log(gamma[2, 2]) + poisson_lpmf(y[i, t] | e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * log(y[i, t-1] + 1) + B[2] * (z_it[i, t] + 1)));
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
    return 0;
  }
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  array[ndept, time] int y;           // data matrix
 // vector[time] r;                  // Trend component
 // vector[time] s;                  // Seasonal component
  //vector[ndept] u;                 // Spatial component  
  simplex[nstate] init_density;       // initial distribution of the Markov chain
  matrix[ndept, time] e_it;           // initial Susceptibles
  matrix[ndept, ndept] R;             // Structure / Precision matrix (IGMRF1/IGMRF2)
  int<lower=0, upper=6> Model;        // Model's functional form
  array[ndept, ndept] int adjacencyMatrix;     // Adjacency matrix
}

parameters {
  real<lower=0, upper=1> G12;            // transition to hyperendemic
  real<lower=0, upper=1> G21;            // transition to endemic
  real<lower=0> kappa_u;                 // spatial precision parameter
  real<lower=0> kappa_r;                 // trend precision parameter
  real<lower=0> kappa_s;                 // seasonal precision parameter
  vector[ndept-1] u;                     // Spatial components
  vector[time] r;                        // Trend components
  vector[time] s;                        // Seasonal components
  vector[2] B;                           // autoregressive coefficients
}

transformed parameters {
    real sumC = sum(u[1:(ndept-1)]);
    vector[ndept] uconstrained;
    uconstrained = append_row(u[1:(ndept-1)], -sumC);
}

model {
  // Priors
  G12 ~ beta(1, 1);
  G21 ~ beta(1, 1);
  kappa_u ~ gamma(1, 0.01);
  kappa_r ~ gamma(1, 0.0001);  
  kappa_s ~ gamma(1, 0.0001);
  
  uconstrained ~ IGMRF1(kappa_u, R);
  r ~ randomwalk2(kappa_r);
  s ~ seasonalComp(kappa_s);

  // Likelihood
  target += Stan_Loglikelihood(y, r, s, uconstrained, G(G12, G21), init_density, e_it, B, Model, adjacencyMatrix);
}
