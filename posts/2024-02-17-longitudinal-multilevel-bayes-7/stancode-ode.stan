//
// Stan code for fitting an ODE model to time-series data
//

// this code block defines the ODE model
    // I'm removing the last letter from each variable name
    // just to avoid potential conflicts with bulit-in names of Stan
    // such as beta for a beta distribution
    // it might not matter, but I wanted to be safe
functions {
  vector odemod(real t,
             vector y,
             real alph, 
             real bet, 
             real gamm, 
             real et) {
    vector[3] dydt;
    dydt[1] = - exp(bet) * y[1] * y[3];
    dydt[2] = exp(bet) * y[1] * y[3] - exp(gamm) * y[2];
    dydt[3] = exp(alph) * y[2] - exp(et) * y[3];
    return dydt;
  }
}

data{
   int<lower = 1> Ntot; //number of observations for each individual
   int<lower = 1> Nind; //number of individuals
   array[Nind] int Nobs; //number of observations for each individual
   array[Ntot] real outcome; //virus load
   array[Ntot] real time; // times at which virus load is measured
   real tstart; //starting time for model
   array[Ntot] int id;  //vector of person IDs to keep track which data points belong to whom
   array[Nind] real dose_adj; //dose after adjustment, 1 value per individual
   array[Nind] int dose_level; //dose level for each individual, needed to index V0 starting values
   //everything below are variables that contain values for prior distributions
   real mu_a_mu; 
   real mu_b_mu;
   real mu_g_mu;
   real mu_e_mu;
   real mu_a_sd;
   real mu_b_sd;
   real mu_g_sd;
   real mu_e_sd;
}

// specifying where in the vector each individual starts and stops
transformed data {
  array[Nind] int start;
  array[Nind] int stop;
  start[1] = 1;
  stop[1] = Nobs[1];
  for(i in 2:Nind) {
    start[i] = start[i - 1] + Nobs[i - 1];
    stop[i] = stop[i - 1] + Nobs[i];
  }
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
    array[Nind] real a0;
    array[Nind] real b0;
    array[Nind] real g0;
    array[Nind] real e0;
    // starting value of virus for individuals in each dose group
    // is being estimated
    array[3] real V0;
}

// Generated/intermediate parameters
transformed parameters{

    // time series for all ODE model variables
    // in our example, each individual has the same number of data points
    // and all data points are collected at the same time
    // if this is not the case, one would need to recode this
    // see e.g. here for an example where the data is more complex
    // https://blog.djnavarro.net/posts/2023-06-10_pop-pk-models/
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
    // predicted virus load from model
    vector[Ntot] virus_pred; 
    array[Ntot] vector[3] y_all;
    vector[3] ystart;

    // loop over all individuals
    for ( i in 1:Nind ) {

      // compute main model parameters
      alph[i] = a0[i] + a1 * dose_adj[i];
      bet[i] = b0[i] + b1 * dose_adj[i];
      gamm[i] = g0[i] + g1 * dose_adj[i];
      et[i] = e0[i] + e1 * dose_adj[i];
     
      // starting value for virus depends on dose 
      ystart = [1e8,0,exp(V0[dose_level[i]])]';
     
     // run ODE for each individual
      y_all[start[i]:stop[i]] = ode_bdf(
        odemod,      // name of ode function
        ystart,      // initial state
        tstart,      // initial time
        time[start[i]:stop[i]],  // observation times - here same for everyone
        alph[i], bet[i], gamm[i], et[i] // model parameters
        );
    
      for (j in 1:Nobs[i]) {
          virus_pred[start[i] + j - 1] = y_all[start[i] + j - 1, 3];
      }
    } // end loop over each individual    
} // end transformed parameters block


model{

        // residual population variation
    sigma ~ exponential( 1 ); 
    // variance of priors
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    sigma_e ~ exponential( 1 );
    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 0 , 0.2); 
    b1 ~ normal( 0 , 0.2);
    g1 ~ normal( 0, 0.2 );
    e1 ~ normal( 0 , 0.2 );
    // hyper-priors to allow for adaptive pooling among individuals 
    // values for the distributions are passed into the Stan code as part of the data
    mu_a ~ normal( mu_a_mu , mu_a_sd );
    mu_b ~ normal( mu_b_mu , mu_b_sd );
    mu_g ~ normal( mu_g_mu , mu_g_sd );
    mu_e ~ normal( mu_e_mu , mu_e_sd );
    // individual variation of each ODE model parameter
    a0 ~ normal( mu_a , sigma_a );
    b0 ~ normal( mu_b , sigma_b );
    g0 ~ normal( mu_g , sigma_g );
    e0 ~ normal( mu_e , sigma_e );
    // prior for starting value 
    V0 ~ normal(5,1);

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
    real<lower=0> sigma_prior;
    real<lower=0> sigma_a_prior;
    real<lower=0> sigma_b_prior;
    real<lower=0> sigma_g_prior;
    real<lower=0> sigma_e_prior;
    real a1_prior;
    real b1_prior;
    real g1_prior;
    real e1_prior;
    real mu_a_prior;
    real mu_b_prior;
    real mu_g_prior;
    real mu_e_prior;
    real a0_prior;     // same prior for each individual so only specify one
    real b0_prior;     
    real g0_prior;     
    real e0_prior;     
    real V0_prior;     
    
    
    // this is so one can plot priors and compare with posterior later   
    // simulate the priors
    sigma_prior = exponential_rng( 1 );
    sigma_a_prior =  exponential_rng(  1 );
    sigma_b_prior = exponential_rng(  1 );
    sigma_g_prior = exponential_rng(  1 );
    sigma_e_prior = exponential_rng(  1 );
    a1_prior = normal_rng( 0 , 0.2);
    b1_prior = normal_rng( 0 , 0.2);
    g1_prior = normal_rng( 0 , 0.2);
    e1_prior = normal_rng( 0 , 0.2);
    mu_a_prior = normal_rng( mu_a_mu , mu_a_sd );
    mu_b_prior = normal_rng( mu_b_mu , mu_b_sd);
    mu_g_prior = normal_rng( mu_g_mu , mu_g_sd);
    mu_e_prior = normal_rng( mu_e_mu , mu_e_sd);
  
    a0_prior = normal_rng(mu_a, sigma_a);
    b0_prior = normal_rng(mu_b, sigma_b);
    g0_prior = normal_rng(mu_g, sigma_g);
    e0_prior = normal_rng(mu_e, sigma_e);

    V0_prior = normal_rng( 5 , 1);


  // compute log-likelihood and predictions
    for(i in 1:Ntot)
    {
      log_lik[i] = normal_lpdf(outcome[i] | virus_pred[i], sigma);
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }
} //end generated quantities block 
