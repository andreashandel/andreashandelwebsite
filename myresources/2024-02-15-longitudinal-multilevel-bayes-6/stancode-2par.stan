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
    mu_a ~ normal( 2 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    sigma_a ~ cauchy( 0 , 1 );
    sigma_b ~ cauchy( 0 , 1 );
  
    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
  
    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 0.1 , 0.1); 
    b1 ~ normal( -0.1 , 0.1);
  
    // residual population variation
    sigma ~ cauchy( 0 , 1 ); 
    
    outcome ~ normal( virus_pred, sigma );
} // close model block

// generated quantities {
//   
//   
// }
