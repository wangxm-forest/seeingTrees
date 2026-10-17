##Copied from gpme repo##


functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1, array[] real x1,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =    gp_exp_quad_cov(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 =  gp_exp_quad_cov(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =    gp_exp_quad_cov(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  vector[N_obs] y_obs;
  
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}

parameters {
  real<lower=0> rho;
  real<lower=0> gamma;
  real<lower=0> sigma;
}

model {
  matrix[N_obs, N_obs] cov =  gp_exp_quad_cov(x_obs, gamma, rho)
                             + diag_matrix(rep_vector(square(sigma), N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(cov);
  
  gamma ~ normal(0,1);
  rho ~ normal(5, 2);
  sigma ~ normal(0,1);

  y_obs ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}

generated quantities {
  vector[N_pred] f_pred = gp_pred_rng(x_pred, y_obs, x_obs, gamma, rho, sigma, 1e-10);
  array[N_pred] real y_pred = normal_rng(f_pred, sigma);
}
