## //
## // Stan code for fitting an ODE model to time-series data
## //

functions {
  vector odemod(real t,
             vector y,
             real alpha, real beta, real gamma, real eta) {
    vector[3] dydt;
    dydt[1] = - beta * y[1] * y[3];
    dydt[2] = beta * y[1] * y[3] - gamma * y[2];
    dydt[3] = alpha * y[2] - eta * y[3];
    return dydt;
  }
}

data{
   int<lower = 1> Nobs; //number of observations
   int<lower = 1> Nind; //number of individuals
   vector[Nobs] outcome;
   vector[Nobs] time;
   vector[Nobs] dose_adj;
   array[Nobs] int id;
}

parameters{
     vector[Nind] a0;
     vector[Nind] b0;
     real mu_a;
     real mu_b;
     real<lower=0> sigma_a;
     real<lower=0> sigma_b;
     real a1;
     real b1;
     real<lower=0> sigma;
}

model{
    
    vector[Nobs] mu;
    vector[Nobs] alpha;
    vector[Nobs] beta;
    
    sigma ~ cauchy( 0 , 1 );
    b1 ~ normal( -0.3 , 1 );
    a1 ~ normal( 0.3 , 1 );
    sigma_b ~ cauchy( 0 , 1 );
    sigma_a ~ cauchy( 0 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    mu_a ~ normal( 2 , 1 );
    b0 ~ normal( mu_b , sigma_b );
    a0 ~ normal( mu_a , sigma_a );
    
    for ( i in 1:Nobs ) {
        beta[i] = b0[id[i]] + b1 * dose_adj[i];
    }
    
    for ( i in 1:Nobs ) {
        alpha[i] = a0[id[i]] + a1 * dose_adj[i];
    }
    
    // deterministic model, given by ODE
    for ( i in 1:Nobs ) {
        mu[i] = exp(alpha[i]) * log(time[i]) - exp(beta[i]) * time[i];
    }
    
    outcome ~ normal( mu , sigma );
}
generated quantities{
    vector[Nobs] log_lik;
     vector[Nobs] mu;
     vector[Nobs] alpha;
     vector[Nobs] beta;
    for ( i in 1:Nobs ) {
        beta[i] = b0[id[i]] + b1 * dose_adj[i];
    }
    for ( i in 1:Nobs ) {
        alpha[i] = a0[id[i]] + a1 * dose_adj[i];
    }
    for ( i in 1:Nobs ) {
        mu[i] = exp(alpha[i]) * log(time[i]) - exp(beta[i]) * time[i];
    }
    for ( i in 1:Nobs ) log_lik[i] = normal_lpdf( outcome[i] | mu[i] , sigma );
}