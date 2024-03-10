//
// Stan code for simulating two ODE models 
//

// the functions code block defines the 2 ODE models
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
    dydt[1] = - bet * exp(y[3]);
    dydt[2] = bet * exp(y[1] + y[3] - y[2]) - gamm;
    dydt[3] = alph * exp(y[2] - y[3]) - et;
    return dydt;
  }
}

data{
   int<lower = 1> Ntot; //number of observations for each individual
   int<lower = 1> Nind; //number of individuals
   int<lower = 1> Ndose; //number of dose levels
   array[Nind] int Nobs; //number of observations for each individual
   array[Nind] int start; //start index for observations for each individual
   array[Nind] int stop; //stop index for observations for each individual
   array[Ntot] real outcome; //virus load
   array[Ntot] real time; // times at which virus load is measured
   real tstart; //starting time for model
   array[Ntot] int id;  //vector of person IDs to keep track which data points belong to whom
   array[Ntot] real dose_adj; //dose after adjustment, 1 value per individual
   array[Ntot] int dose_level; //dose level for each individual, needed to index V0 starting values
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

// everything that is estimated goes here
parameters{
     // individual-level  parameters
    // right now just the non-exponentiated form of the main ones
    vector[Nind] a0;
    vector[Nind] b0;
    vector[Nind] g0;
    vector[Nind] e0;
    // starting value of virus for individuals in each dose group
    // is being estimated
    // vector[Ndose] V0;
    // population variance
    real<lower=0> sigma;

    }


transformed parameters{

    // main model parameters
    // this is coded such that each individual can have their own value
    vector[Nind] alph;
    vector[Nind] bet;
    vector[Nind] gamm;
    vector[Nind] et;

    // time series for all 3 ODE model variables
    array[Ntot] vector[3] y_all;

    // predicted virus load from model
    vector[Ntot] virus_pred; 


    // starting conditions for ODE model
    vector[3] ystart;
    
      // starting value for virus depends on dose 
      // we are fitting/running model with variables on a log scale
      ystart =  [ 1e12, 1e-15, 1e-3]';

    // loop over all individuals
    for ( i in 1:Nind ) {


      // compute main model parameters
      // here just exponentiated version of estimated parameters
      alph[i] = exp(a0[i]) ;
      bet[i] =  exp(b0[i]) ;
      gamm[i] =  exp(g0[i]) ;
      et[i] =  exp(e0[i]) ;
     
 
     // repeat ODE model
      y_all[start[i]:stop[i]] = ode_rk45(
        odemod,      // name of ode function
        log(ystart),      // initial state
        tstart,      // initial time
        time[start[i]:stop[i]],  // observation times - here same for everyone
        alph[i], // model parameters - exponentiated to enforce positivity 
        bet[i], 
        gamm[i], 
        et[i] 
        );
      virus_pred[start[i]:stop[i]] = to_vector(y_all[start[i]:stop[i], 3]);

    } // end loop over each individual    
} // end transformed parameters block


// anything that is probabilistic goes here
model{

        // residual population variation
    sigma ~ cauchy(0, 1 ); 
    // average dose-dependence of each ODE model parameter
    a0 ~ normal( a0_mu , a0_sd); 
    b0 ~ normal( b0_mu , b0_sd);
    g0 ~ normal( g0_mu,  g0_sd);
    e0 ~ normal( e0_mu , e0_sd);

    // distribution of outcome (virus load)
    // all computations to get the time-series trajectory for the outcome are done  
    // inside the transformed parameters block
    outcome ~ normal( virus_pred , sigma );

  }

// for model diagnostics and exploration
generated quantities {

    // define quantities that are computed in this block
    vector[Ntot] ypred;

  // compute log-likelihood and predictions
    for(i in 1:Ntot)
    {
      ypred[i] = normal_rng(virus_pred[i], sigma);
    }


} //end generated quantities block 

