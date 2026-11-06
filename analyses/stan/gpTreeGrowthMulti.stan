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
  int<lower=1> N_tree; //Number of individual trees
  int<lower=1> N_obs;
  int<lower=1> N_all_year_obs;
  
  vector[N_obs] temp_obs;
  
  array[N_tree] int<lower=1> start;
  array[N_tree] int<lower=1> end;
  
  array[N_tree] int<lower=1, upper=N_obs> N_obs_tree; # Number of Observations per Tree
  array[N_all_year_obs] int all_year_obs; # List of all Years Used
  
  # Year and rw observations
  real year_obs[N_obs]; 
  vector[N_obs] rw_obs;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> rho;
  real<lower=0> gamma;
  real<lower=0> sigma;
}

transformed parameters{
  vector[N_obs] y_clim;
  for (n in 1:N_obs){
    y_clim[n] = alpha + beta * temp_obs[n];
  }
}

model {
  alpha ~ normal(0.5, 0.15);
  beta ~ normal(0, 0.5);
  gamma ~ normal(0, 1);
  rho ~ normal(100, 50);
  sigma ~ normal(0, 1);
  
# Compute the covariance matrix across all years
  matrix[N_all_year_obs, N_all_year_obs] cov = gp_exp_quad_cov(all_year_obs, gamma, rho)
                                    + diag_matrix(rep_vector(square(sigma), N_all_year_obs));
  matrix[N_all_year_obs, N_all_year_obs] L_cov = cholesky_decompose(cov);
  
  
 # Subset L_cov for each tree
  for (i in 1:N_tree){
    matrix[N_obs_tree[i], N_obs_tree[i]] L_cov_tree = block(L_cov, 1, 1, N_obs_tree[i], N_obs_tree[i]);
    rw_obs[start[i]:end[i]] ~ multi_normal_cholesky(y_clim[start[i]:end[i]], L_cov_tree);
  }
}

generated quantities {
  vector[N_obs] f;
  array[N_obs] real rw_pred;
  
  for (i in 1:N_tree){
    f[start[i]:end[i]] = gp_pred_rng(year_obs[start[i]:end[i]], 
                                    rw_obs[start[i]:end[i]] - y_clim[start[i]:end[i]], 
                                    year_obs[start[i]:end[i]],
                                    gamma, rho, sigma, 1e-10);
  
    rw_pred[start[i]:end[i]] = normal_rng(f[start[i]:end[i]], sigma);
  }
}
