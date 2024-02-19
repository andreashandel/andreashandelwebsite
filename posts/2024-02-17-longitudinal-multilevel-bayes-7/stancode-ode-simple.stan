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
    // population-level dose dependence parameters
    real a1;
    real b1;
    real g1;
    real e1;
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
      alph[i] =  a1 * (dose_adj[i] + 1);
      bet[i] =  b1 * (dose_adj[i] + 1);
      gamm[i] =  g1 * (dose_adj[i] + 1);
      et[i] =  e1 * (dose_adj[i] + 1);
     
      // starting value for virus depends on dose 
      ystart = [1e8,0,exp(V0[dose_level[i]])]';
     
     // run ODE for each individual
      y_all[start[i]:stop[i]] = ode_rk45(
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
    // average dose-dependence of each ODE model parameter
    a1 ~ normal( 10 , 1); 
    b1 ~ normal( -23 , 1);
    g1 ~ normal( 1,  1);
    e1 ~ normal( 1 , 1);
    // prior for starting value 
    V0 ~ normal(5,1);

    // distribution of outcome (virus load)
    // all computations to get the time-series trajectory for the outcome are done  
    // inside the transformed parameters block
    outcome ~ normal( virus_pred , sigma );

}
