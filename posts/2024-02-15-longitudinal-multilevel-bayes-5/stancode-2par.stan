//
// Stan code for fitting a hierarchical model to time-series data
//

data{
   int<lower = 1> Nobs; //number of observations
   int<lower = 1> Nind; //number of individuals
   vector[Nobs] outcome; //virus load
   vector[Nobs] time; // times at which virus load is measured
   vector[Nobs] dose_adj; //dose after adjustment, only differs by individual but is repeated
   array[Nobs] int id;  //vector of person IDs to keep track which data points belong to whom
}

parameters{
     //  population variance
     real<lower=0> sigma;
     // population-level dose dependence parameters
     real a1;
     real b1;
     // individual level variation parameters
     vector[Nind] a0;
     vector[Nind] b0;
     // hyper-parameters to implement adaptive pooling
     real mu_a;
     real mu_b;
     real<lower=0> sigma_a;
     real<lower=0> sigma_b;
} // close parameters block


// Generated/intermediate parameters
transformed parameters{
    // predicted virus load from model
    vector[Nobs] virus_pred; 
    // main model parameters
    // each individual has their potentially own value 
    vector[Nind] alpha;
    vector[Nind] beta;
 
    // compute main model parameters
    // note that dos_adj contains a dose for each observation, 
    // but it's repeated after the first Nind entries, so we don't need to loop over all obersvations
    for ( i in 1:Nind ) {
        alpha[i] = a0[i] + a1 * dose_adj[i];
        beta[i] = b0[i] + b1 * dose_adj[i];
    }
    // loop over all observations
    // since paramaters are saved in vectors of length corresponding to number of individuals
    // we need to index with that extra id[i] notation
    for (i in 1:Nobs)
    {
      virus_pred[i] = exp(alpha[id[i]]) * log(time[i]) - exp(beta[id[i]]) * time[i] ;
    }
} // close transformed parameters block



model{

    // hyper-priors to allow for adaptive pooling among individuals 
    mu_a ~ normal( 3 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    sigma_a ~ exponential(  1 );
    sigma_b ~ exponential(  1 );
  
    // individual variation of each model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
  
    // average dose-dependence of each model parameter
    a1 ~ normal( 0.1 , 0.1); 
    b1 ~ normal( -0.1 , 0.1);
  
    // residual population variation
    sigma ~ exponential(  1 ); 
    
    outcome ~ normal( virus_pred, sigma );
} // close model block

// this code block is not needed for the model fitting
// but it's useful to compute quantities that one can use later
// for model diagnostics and exploration
generated quantities {
    // define quantities that are computed in this block
    vector[Nobs] ypred;
    vector[Nobs] log_lik;
    real<lower=0> sigma_prior;
    real<lower=0> sigma_a_prior;
    real<lower=0> sigma_b_prior;
    real a1_prior;
    real b1_prior;
    real mu_a_prior;
    real mu_b_prior;
    real a0_prior;     // same prior for each individual so only specify one
    real b0_prior;     
    
    
    // this is so one can plot priors and compare with posterior later   
    // simulate the priors
    sigma_prior = exponential_rng(  1 );
    sigma_a_prior =  exponential_rng( 1 );
    sigma_b_prior = exponential_rng(  1 );
    a1_prior = normal_rng( 0.1 , 0.1);
    b1_prior = normal_rng( -0.1 , 0.1);
    mu_a_prior = normal_rng( 3 , 1 );
    mu_b_prior = normal_rng( 0.5 , 1 );
    a0_prior = normal_rng(mu_a, sigma_a);
    b0_prior = normal_rng(mu_b, sigma_b);
  
  // compute log-likelihood and predictions
    for(i in 1:Nobs)
    {
      log_lik[i] = normal_lpdf(outcome[i] | virus_pred[i], sigma);
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }
} //end generated quantities block 
