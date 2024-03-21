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
   real a0_mu; 
   real b0_mu;
   real g0_mu;
   real e0_mu;
   real a0_sd;
   real b0_sd;
   real g0_sd;
   real e0_sd;
}

parameters{
    // population variance
    real<lower=0> sigma;
    // individual-level  parameters
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
    for ( i in 1:Nind ) {
        alph[i] = exp(30*a0[i]);
        bet[i] = exp(b0[i]);
        gamm[i] = exp(g0[i]);
        et[i] = exp(e0[i]);
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

    // individual variation of each ODE model parameter
    a0 ~ normal( a0_mu , a0_sd );
    b0 ~ normal( b0_mu , b0_sd );
    g0 ~ normal( g0_mu , g0_sd );
    e0 ~ normal( e0_mu , e0_sd );

    // residual population variation
    sigma ~ cauchy(0, 1); 

    // distribution of outcome (virus load)
    // all computations to get the time-series trajectory for the outcome are done  
    // inside the transformed parameters block
    outcome ~ normal( virus_pred , sigma );
}

// for model diagnostics and exploration
generated quantities {
    // define quantities that are computed in this block
    vector[Ntot] ypred;
    vector[Ntot] log_lik;
    //real<lower=0> sigma_prior;
    real sigma_prior;
    real a0_prior;
    real b0_prior;
    real g0_prior;
    real e0_prior;
        
    // this is so one can plot priors and compare with posterior later   
    // simulate the priors
    sigma_prior = abs(cauchy_rng(0, 1));
    a0_prior = normal_rng( a0_mu , a0_sd);
    b0_prior = normal_rng( b0_mu , b0_sd);
    g0_prior = normal_rng( g0_mu , g0_sd);
    e0_prior = normal_rng( e0_mu , e0_sd);


  // compute log-likelihood and predictions
    for(i in 1:Ntot)
    {
      log_lik[i] = normal_lpdf(outcome[i] | virus_pred[i], sigma);
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }
} //end generated quantities block 
