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
   vector[Nobs] outcome; //virus load
   vector[Nobs] time;
   vector[Nobs] dose_adj; //dose after adjustment
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
    
    // predicted virus load from model
    vector[Nobs] virus_pred; 
    // computed quantities for each ODE model parameter
    vector[Nobs] alpha;
    vector[Nobs] beta;
    vector[Nobs] gamma;
    vector[Nobs] eta;
    
    // residual population variation
    sigma ~ cauchy( 0 , 1 ); 
    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 0.3 , 1 ); 
    b1 ~ normal( -0.3 , 1 );
    g1 ~ normal( -0.3 , 1 );
    e1 ~ normal( -0.3 , 1 );
    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
    g0 ~ normal( mu_g , sigma_g );
    e0 ~ normal( mu_e , sigma_e );
    // hyper-priors to allow for adaptive pooling among individuals 
    mu_a ~ normal( 2 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    mu_g ~ normal( 0.5 , 1 );
    mu_e ~ normal( 0.5 , 1 );
    sigma_a ~ cauchy( 0 , 1 );
    sigma_b ~ cauchy( 0 , 1 );
    sigma_g ~ cauchy( 0 , 1 );
    sigma_e ~ cauchy( 0 , 1 );
    
    
    for ( i in 1:Nobs ) {
        alpha[i] = a0[id[i]] + a1 * dose_adj[i];
        beta[i] = b0[id[i]] + b1 * dose_adj[i];
        gamma[i] = g0[id[i]] + g1 * dose_adj[i];
        eta[i] = e0[id[i]] + e1 * dose_adj[i];
        
    }
    
    for ( i in 1:Nind ) {
      virus_pred = ode_rk45(
      amount_change,            // ode function
      initial_amount[i],        // initial state
      initial_time,             // initial time
      t_obs[start[i]:stop[i]],  // observation times
      KA[i],                    // absorption rate
      CL[i],                    // clearance
      V[i]                      // volume
    );
    }
    
    outcome ~ normal( log(virus_pred) , sigma );
}

