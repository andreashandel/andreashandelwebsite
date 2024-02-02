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
    // since paramaters are saved in vectors of length corresponding to number of individuals
    // we need to index with that extra id[i] notation
    for (i in 1:Nobs)
    {
      virus_pred[i] = eta[id[i]] + exp(alpha[id[i]]) * log(time[i]) - exp(beta[id[i]]) * (time[i] ^ exp(gamma[id[i]]));
    }

}



model{

    // hyper-priors to allow for adaptive pooling among individuals 
    mu_a ~ normal( 2 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    mu_g ~ normal( 0 , 0.01 );
    mu_e ~ normal( 1 , 1 );
    sigma_a ~ cauchy( 0 , 1 );
    sigma_b ~ cauchy( 0 , 1 );
    sigma_g ~ cauchy( 0 , 1 );
    sigma_e ~ cauchy( 0 , 1 );

    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
    g0 ~ normal( mu_g , sigma_g );
    e0 ~ normal( mu_e , sigma_e );

    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 0.1 , 0.1); 
    b1 ~ normal( -0.1 , 0.1);
    g1 ~ normal( 0 , 0.01 );
    e1 ~ normal( 0 , 0.1 );

    // residual population variation
    sigma ~ cauchy( 0 , 1 ); 
    
    outcome ~ normal( log(virus_pred) , sigma );
}

generated quantities {
  
  
}
