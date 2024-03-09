//
// Stan code for simulating two ODE models 
//

// the functions code block defines the 2 ODE models
    // I'm removing the last letter from each variable name
    // just to avoid potential conflicts with bulit-in names of Stan
    // such as beta for a beta distribution
    // it might not matter, but I wanted to be safe
functions {
  vector odemod1(real t,
             vector y,
             real alph, 
             real bet, 
             real gamm, 
             real et) {
    vector[3] dydt;
    dydt[1] = - bet * y[1] * y[3];
    dydt[2] = bet * y[1] * y[3] - gamm * y[2];
    dydt[3] = alph * y[2] - et * y[3];
    return dydt;
  }

  vector odemod2(real t,
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
   real V0_mu;
   real V0_sd;
}


parameters{
    }


model{
  }

// for model diagnostics and exploration
generated quantities {
    // main model parameters
    // this is coded such that each individual can have their own value
    vector[Nind] alph;
    vector[Nind] bet;
    vector[Nind] gamm;
    vector[Nind] et;

    // individual-level  parameters
    // right now just the non-exponentiated form of the main ones
    vector[Nind] a0;
    vector[Nind] b0;
    vector[Nind] g0;
    vector[Nind] e0;

  
    // starting value of virus for individuals in each dose group
    // is being estimated
    vector[Ndose] V0;
    // population variance
    real<lower=0> sigma;

    // time series for all 3 ODE model variables
    array[Ntot] vector[3] y_all1;
    array[Ntot] vector[3] y_all2;

    // predicted virus load from model
    vector[Ntot] virus_pred1; 
    vector[Ntot] virus_pred2; 

    // define quantities that are computed in this block
    vector[Ntot] ypred1;
    vector[Ntot] ypred2;

    // starting conditions for ODE model
    vector[3] ystart;

    // compute fitness of virus
    real R0;

    // sample residual variation
    sigma = abs(cauchy_rng(0, 1 ));
    
    // loop over dose for starting value
    for (n in 1:3)
    {
        V0[n] = exp(normal_rng(V0_mu, V0_sd));
    }

    // loop over all individuals
    for ( i in 1:Nind ) {


    a0[i] = normal_rng( a0_mu , a0_sd);
    b0[i] = normal_rng( b0_mu , b0_sd);
    g0[i] = normal_rng( g0_mu , g0_sd);
    e0[i] = normal_rng( e0_mu , e0_sd);


      // compute main model parameters
      // here just exponentiated version of estimated parameters
      alph[i] = exp(a0[i]) ;
      bet[i] =  exp(b0[i]) ;
      gamm[i] =  exp(g0[i]) ;
      et[i] =  exp(e0[i]) ;
     
      // starting value for virus depends on dose 
      // we are fitting/running model with variables on a log scale
      ystart =  [ 1e8, 1, V0[dose_level[start[i]]] ]';
 
     // run ODE for each individual
      y_all1[start[i]:stop[i]] = ode_bdf(
        odemod1,      // name of ode function
        ystart,      // initial state
        tstart,      // initial time
        time[start[i]:stop[i]],  // observation times - here same for everyone
        alph[i], // model parameters - exponentiated to enforce positivity 
        bet[i], 
        gamm[i], 
        et[i] 
        );
      virus_pred1[start[i]:stop[i]] = log(abs(to_vector(y_all1[start[i]:stop[i], 3])));

     // repeat for second ODE model
      y_all2[start[i]:stop[i]] = ode_bdf(
        odemod2,      // name of ode function
        log(ystart),      // initial state
        tstart,      // initial time
        time[start[i]:stop[i]],  // observation times - here same for everyone
        alph[i], // model parameters - exponentiated to enforce positivity 
        bet[i], 
        gamm[i], 
        et[i] 
        );
      virus_pred2[start[i]:stop[i]] = to_vector(y_all2[start[i]:stop[i], 3]);

    } // end loop over each individual    

  // compute log-likelihood and predictions
    for(i in 1:Ntot)
    {
      ypred1[i] = normal_rng(virus_pred1[i], sigma);
      ypred2[i] = normal_rng(virus_pred2[i], sigma);
    }


} //end generated quantities block 

