//
// Stan code for fitting a hierarchical model to time-series data
//

data{
   int<lower = 1> Ntot; //number of observations
   int<lower = 1> Nind; //number of individuals
   vector[Ntot] outcome; //virus load
   vector[Ntot] time; // times at which virus load is measured
   vector[Nind] dose_adj; //dose after adjustment, 1 value per individual
   array[Ntot] int id;  //vector of person IDs to keep track which data points belong to whom
   //everything below are variables that contain values for prior distributions
   real mu_a_mu; 
   real mu_b_mu;
   real mu_g_mu;
   real mu_e_mu;
   real mu_a_sd;
   real mu_b_sd;
   real mu_g_sd;
   real mu_e_sd;
   real a1_mu; 
   real b1_mu;
   real g1_mu;
   real e1_mu;
   real a1_sd;
   real b1_sd;
   real g1_sd;
   real e1_sd;

}

parameters{
    // population variance
    real<lower=0> sigma;
    // variance of priors
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> sigma_g;
    real<lower=0> sigma_e;
    // population-level dose dependence parameters
    real a1;
    real b1;
    real g1;
    real e1;
    // hyper-parameters to implement adaptive pooling
    real mu_a;
    real mu_b;
    real mu_g;
    real mu_e;
    // individual level variation parameters
    vector[Nind] a0;
    vector[Nind] b0;
    vector[Nind] g0;
    vector[Nind] e0;
}


// Generated/intermediate parameters
transformed parameters{

    // predicted virus load from model
    vector[Ntot] virus_pred; 
    // main model parameters
    // each individual has their potentially own value 
    // I'm removing the last letter from each variable name
    // just to avoid potential conflicts with bulit-in names of Stan
    // such as beta for a beta distribution
    // it might not matter, but I wanted to be safe
    vector[Nind] alph;
    vector[Nind] bet;
    vector[Nind] gamm;
    vector[Nind] et;
    // intermediate variables for denominators
    // just to make equation easier to understand
    real d1; 
    real d2; 


    // compute main model parameters
    // multiply a by 30 to get better scaling
    for ( i in 1:Nind ) {
        alph[i] = exp( 30*( a0[i] + a1 * dose_adj[i]) );
        bet[i]  = exp( b0[i] + b1 * dose_adj[i] );
        gamm[i] = exp( g0[i] + g1 * dose_adj[i] );
        et[i]   = exp( e0[i] + e1 * dose_adj[i] );
        
    }
    // loop over all observations
    // since parameters are saved in vectors of length corresponding to number of individuals
    // we need to index with that extra id[i] notation
    for (i in 1:Ntot)
    {
      d1 = exp( -bet[id[i]] * (time[i] - gamm[id[i]]) );
      d2 = exp(   et[id[i]] * (time[i] - gamm[id[i]]) ); 
      virus_pred[i] = log(  2*alph[id[i]] / ( d1 + d2) );
     }

} // end transformed parameters block



model{

    // hyper-priors to allow for adaptive pooling among individuals 
    // values for the distributions are passed into the Stan code as part of the data
    mu_a ~ normal( mu_a_mu , mu_a_sd );
    mu_b ~ normal( mu_b_mu , mu_b_sd );
    mu_g ~ normal( mu_g_mu , mu_g_sd );
    mu_e ~ normal( mu_e_mu , mu_e_sd );
    // variance of priors
    sigma_a ~ exponential(1);
    sigma_b ~ exponential(1);
    sigma_g ~ exponential(1);
    sigma_e ~ exponential(1);
    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
    g0 ~ normal( mu_g , sigma_g );
    e0 ~ normal( mu_e , sigma_e );
    // average dose-dependence of each ODE model parameter
    // values for the distributions are passed into the Stan code as part of the data
    a1 ~ normal( a1_mu , a1_sd); 
    b1 ~ normal( b1_mu , b1_sd);
    g1 ~ normal( g1_mu , g1_sd);
    e1 ~ normal( e1_mu , e1_sd);
    // residual population variation
    sigma ~ exponential(1); 

    // distribution of outcome (virus load)
    // all computations to get the time-series trajectory for the outcome are done  
    // inside the transformed parameters block
    outcome ~ normal( virus_pred , sigma );
}

// for model diagnostics and exploration
generated quantities {
    // define quantities that are computed in this block
    // for the individual level parameters, every individual
    // has the same prior so we only specify one
    real mu_a_prior;
    real mu_b_prior;
    real mu_g_prior;
    real mu_e_prior;
    real<lower=0> sigma_a_prior;
    real<lower=0> sigma_b_prior;
    real<lower=0> sigma_g_prior;
    real<lower=0> sigma_e_prior;
    real a0_prior;  
    real b0_prior;     
    real g0_prior;     
    real e0_prior;     

    real a1_prior;
    real b1_prior;
    real g1_prior;
    real e1_prior;

    real<lower=0> sigma_prior;
    vector[Ntot] log_lik;
    vector[Ntot] ypred;

    
    // this is so one can plot priors and compare with posterior later   
    // simulate the priors
    mu_a_prior = normal_rng( mu_a_mu , mu_a_sd);
    mu_b_prior = normal_rng( mu_b_mu , mu_b_sd);
    mu_g_prior = normal_rng( mu_g_mu , mu_g_sd);
    mu_e_prior = normal_rng( mu_e_mu , mu_e_sd);

    sigma_a_prior =  exponential_rng(  1 );
    sigma_b_prior = exponential_rng(  1 );
    sigma_g_prior = exponential_rng(  1 );
    sigma_e_prior = exponential_rng(  1 );

    a0_prior = normal_rng(mu_a_prior, sigma_a_prior);
    b0_prior = normal_rng(mu_b_prior, sigma_b_prior);
    g0_prior = normal_rng(mu_g_prior, sigma_g_prior);
    e0_prior = normal_rng(mu_e_prior, sigma_e_prior);

    a1_prior = normal_rng( a1_mu , a1_sd);
    b1_prior = normal_rng( b1_mu , b1_sd);
    g1_prior = normal_rng( g1_mu , g1_sd);
    e1_prior = normal_rng( e1_mu , e1_sd);

    sigma_prior = exponential_rng( 1 );

  
  // compute log-likelihood and predictions
    for(i in 1:Ntot)
    {
      log_lik[i] = normal_lpdf(outcome[i] | virus_pred[i], sigma);
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }
} //end generated quantities block 
