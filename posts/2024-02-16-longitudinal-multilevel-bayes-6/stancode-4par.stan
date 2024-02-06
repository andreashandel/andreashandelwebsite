//
// Stan code for fitting a hierarchical model to time-series data
//

data{
   int<lower = 1> Nobs; //number of observations
   int<lower = 1> Nind; //number of individuals
   vector[Nobs] outcome; //virus load
   vector[Nobs] time; // times at which virus load is measured
   vector[Nind] dose_adj; //dose after adjustment, 1 value per individual
   array[Nobs] int id;  //vector of person IDs to keep track which data points belong to whom
   real mu_a_mu; //everything below is values for prior distributions
   real mu_b_mu;
   real mu_g_mu;
   real mu_e_mu;
   real mu_a_sd;
   real mu_b_sd;
   real mu_g_sd;
   real mu_e_sd;
}

parameters{
     //  population variance
     real<lower=0> sigma;
     // population-level dose dependence parameters
     real a1;
     real b1;
     real g1;
     real e1;
     // individual level variation parameters
     vector[Nind] a0;
     vector[Nind] b0;
     vector[Nind] g0;
     vector[Nind] e0;
     // hyper-parameters to implement adaptive pooling
     real mu_a;
     real mu_b;
     real mu_g;
     real mu_e;
     real<lower=0> sigma_a;
     real<lower=0> sigma_b;
     real<lower=0> sigma_g;
     real<lower=0> sigma_e;
}


// Generated/intermediate parameters
transformed parameters{

    // predicted virus load from model
    vector[Nobs] virus_pred; 
    // main model parameters
    // each individual has their potentially own value 
    vector[Nind] alpha;
    vector[Nind] beta;
    vector[Nind] gamma;
    vector[Nind] eta;

    // compute main model parameters
    for ( i in 1:Nind ) {
        alpha[i] = a0[i] + a1 * dose_adj[i];
        beta[i] = b0[i] + b1 * dose_adj[i];
        gamma[i] = g0[i] + g1 * dose_adj[i];
        eta[i] = e0[i] + e1 * dose_adj[i];
        
    }
    // loop over all observations
    // since parameters are saved in vectors of length corresponding to number of individuals
    // we need to index with that extra id[i] notation
    for (i in 1:Nobs)
    {
      virus_pred[i] = log( 2*exp(alpha[id[i]]) / ( exp( -exp(beta[id[i]]) * (exp(gamma[id[i]]) - time[i]) )  +  exp( exp(eta[id[i]]) * (time[i] - exp(gamma[id[i]])) )  ) ) ;
     }

} // end transformed parameters block



model{

    // hyper-priors to allow for adaptive pooling among individuals 
    mu_a ~ normal( mu_a_mu , mu_a_sd );
    mu_b ~ normal( mu_b_mu , mu_b_sd );
    mu_g ~ normal( mu_g_mu , mu_g_sd );
    mu_e ~ normal( mu_e_mu , mu_e_sd );
    sigma_a ~ normal( 0 , 1 );
    sigma_b ~ normal( 0 , 1 );
    sigma_g ~ normal( 0 , 1 );
    sigma_e ~ normal( 0 , 1 );

    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
    g0 ~ normal( mu_g , sigma_g );
    e0 ~ normal( mu_e , sigma_e );

    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 0 , 0.2); 
    b1 ~ normal( 0 , 0.2);
    g1 ~ normal( 0.25, 0.2 );
    e1 ~ normal( 0 , 0.2 );

    // residual population variation
    sigma ~ normal( 0 , 2 ); 
    
    outcome ~ normal( virus_pred , sigma );
}

// for model diagnostics and exploration
generated quantities {
    // define quantities that are computed in this block
    vector[Nobs] ypred;
    vector[Nobs] log_lik;
    real<lower=0> sigma_prior;
    real a1_prior;
    real b1_prior;
    real g1_prior;
    real e1_prior;
    real mu_a_prior;
    real mu_b_prior;
    real mu_g_prior;
    real mu_e_prior;
    real<lower=0> sigma_a_prior;
    real<lower=0> sigma_b_prior;
    real<lower=0> sigma_g_prior;
    real<lower=0> sigma_e_prior;
    real a0_prior;     // same prior for each individual so only specify one
    real b0_prior;     
    real g0_prior;     
    real e0_prior;     
    
    
    // this is so one can plot priors and compare with posterior later   
    // simulate the priors
    sigma_prior = normal_rng( 0 , 2 );
    a1_prior = normal_rng( 0 , 0.2);
    b1_prior = normal_rng( 0 , 0.2);
    g1_prior = normal_rng( 0.25 , 0.2);
    e1_prior = normal_rng( 0 , 0.2);
    mu_a_prior = normal_rng( mu_a_mu , mu_a_sd );
    mu_b_prior = normal_rng( mu_b_mu , mu_b_sd);
    mu_g_prior = normal_rng( mu_g_mu , mu_g_sd);
    mu_e_prior = normal_rng( mu_e_mu , mu_e_sd);
    sigma_a_prior =  normal_rng( 0 , 1 );
    sigma_b_prior = normal_rng( 0 , 1 );
    sigma_g_prior = normal_rng( 0 , 1 );
    sigma_e_prior = normal_rng( 0 , 1 );
  
    a0_prior = normal_rng(mu_a, sigma_a);
    b0_prior = normal_rng(mu_b, sigma_b);
    g0_prior = normal_rng(mu_g, sigma_g);
    e0_prior = normal_rng(mu_e, sigma_e);
  
  // compute log-likelihood and predictions
    for(i in 1:Nobs)
    {
      log_lik[i] = normal_lpdf(outcome[i] | virus_pred[i], sigma);
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }
} //end generated quantities block 
