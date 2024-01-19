functions{

   // Intrinsic GMRF density
    real IGMRF_lpdf(vector u, real kappa_u, matrix R) {
    int n = rows(R);
    return (((n - 1) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, u));
  }
  
  // Intrinsic GMRF density
    real IGMRF_lpdf(vector u, matrix Q) {
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
}