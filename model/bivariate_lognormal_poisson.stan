data{
  int N;
  int<lower=0> home_fouls[N];
  int<lower=0> away_fouls[N];
  int k;
  matrix[N,k] X;
  int<lower=2> N_teams;
  int<lower=1> home_team[N];
  int<lower=1> away_team[N];
}
parameters{
  vector[k] home_beta;
  vector[k] away_beta;
  real<lower=0> home_sd;
  real<lower=0> away_sd;
  real<lower=-1,upper=1> rho;
  matrix[2,N] z;
  real<lower=0> home_team_sd;
  vector[N_teams] home_tilde;
  real<lower=0> away_team_sd;
  vector[N_teams] away_tilde;
}
transformed parameters{
  vector[N_teams] home = home_team_sd*home_tilde;
  vector[N_teams] away = away_team_sd*away_tilde;
  cov_matrix[2] Sigma;
    
  Sigma[1,1] = square(home_sd);
  Sigma[1,2] = rho*home_sd*away_sd;
  Sigma[2,1] = rho*home_sd*away_sd;
  Sigma[2,2] = square(away_sd);
  
}
model{
  matrix[2,2] L_Sigma = cholesky_decompose(Sigma);
  matrix[2,N] Z =  L_Sigma*z;
  
  home_beta[1] ~ normal(2.5, 1);
  home_beta[2:k] ~ normal(0, 1);
  away_beta[1] ~ normal(2.5, 1);
  away_beta[2:k] ~ normal(0, 1);
  
  home_sd ~ normal(0, 2.5);
  away_sd ~ normal(0, 2.5);
  
  0.5*(rho+1) ~ beta(3, 2);
  
  to_vector(z) ~ normal(0, 1);
  
  home_team_sd ~ normal(0, 2.5);
  home_tilde ~ normal(0, 1);
  
  away_team_sd ~ normal(0, 2.5);
  away_tilde ~ normal(0, 1);

  home_fouls ~ poisson_log(X*home_beta + home[home_team] + Z[1,]');
  away_fouls ~ poisson_log(X*away_beta + away[away_team] + Z[2,]');
}
