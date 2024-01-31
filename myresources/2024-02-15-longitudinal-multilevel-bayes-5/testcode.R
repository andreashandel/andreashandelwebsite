library(brms)
library(tidyverse)
library(here)
library(rstan)
library(cmdstanr)

rstan_options(auto_write=TRUE)
#options(mc.cores=parallel::detectCores())

# simulated data for exponential decay: c=A0*exp(-k1*time)
set.seed(1)
time <- c(0,1,2,3,4,5,7,9,11)
y_exp <- 100*exp(-0.2*time)  # c0 =100, kr=0.2, exact simulated data
df_sim <- data.frame(time, y_exp)
df_sim <- df_sim %>% mutate(y_obs=y_exp+rnorm(length(time),0,5)) %>% dplyr::select(-y_exp)# random normal error added to simulated data with mean = 0 and sd = 5
N_t <- length(df_sim$time) # number of datapoints





# the following code uses the new integrators directly in Stan, works also fine

write("
functions {                  // input to the ode_integrator:
  vector fo(real t,          // time
            vector y,        // state
            vector kr        // parameters (only one in this case, c0 is automatically estimated)
                  ) {      
    vector[1] dydt;
    dydt[1] = -kr[1] * y[1];
    return dydt;
  }
}
data {
  int<lower = 1>   N;
  real             t0;
  real<lower = t0> ts[N];      // each element must be greater than t0
  real             y0_obs;     // the observed initial value
  real             y_obs[N];   //the observed (measured) values excluding the initial value
}
parameters {
  real<lower = 0>      sigma; // experimental sigma
  vector<lower = 0>[1] kr;    // rate constant to be estimated
  vector[1]            y0;    // initial condition to be estimated
}
model {
  vector[1] mu[N] = ode_bdf(fo, y0, t0, ts, kr);  // call to the ode_integrator, result is put in mu
  sigma ~ cauchy(0, 5);                           // prior for sigma
  kr[1] ~ normal(0.2, 0.02);                      // prior for rate constant
  y0 ~ normal(100,10);                            // prior on initial condition
  y0_obs ~ normal(y0[1], sigma);
  for (t in 1:N) {
    y_obs[t] ~ normal(mu[t], sigma);
  }
}    

","stan_fo.stan")

# compile the stan model:
fo_model1 <- cmdstan_model("stan_fo.stan")
# the data list for Stan:
fo_data <- list(
  N = N_t - 1,                # the first datapoint should not be at time zero because that point is provided by t0 and y0
  t0 = 0,
  ts = df_sim$time[-1],       # omitting the first time datapoint at time zero because that needs to be given for y0 (next line)
  y0_obs = df_sim$y_obs[1],   # the first datapoint is the initial value y0 
  y_obs = df_sim$y_obs[-1]    # the measured data without the starting point
)

# sampling from the compiled model:
fo_fit <- fo_model1$sample(
  data = fo_data,
  seed = 42L,
  refresh = 0,                
  parallel_chains=4 
)

# So far, so good but when coupled to brms things go wrong:

fo_model2 <- "
  vector fo2(real t,                  //time
             vector y,                // the state
             real k1){                // the parameter next to A0
    
             vector[1] dydt;          // dimension of ODE
             dydt[1] = - k1 * y[1];   // ODE
             return dydt;             // returns a 1-dimensional vector     
  }

// function called from brms: integration of ODEs
vector[1,1]  fo_ode2(t,A0,k1) {
        real<lower=0> k1;
        vector[1] y0;        //one initial value
        y0[1]=A0;            //initial value
        vector [1] y_hat[N]=ode_rk45(fo2, y0,0.0,time,k1);
  return(y[1,1]);
}

"

df_sim_ode <- df_sim[-1,] #remove the first datapoint at t0

fo_formula <- bf(y_obs~fo_ode2(time,A0,k1),
                 A0~1,
                 k1~1,
                 nl=TRUE)

fo_priors <- c(prior(normal(100,10),nlpar=A0),
               prior(normal(0.2,0.1), nlpar=k1), 
               prior(cauchy(0,10), class=sigma)
)

fo_result2 <- brm(data=df_sim_ode,family = gaussian,formula=fo_formula,prior=fo_priors, inits=0,iter=4000,chains=4, cores=4, stanvars = stanvar(scode=fo_model2, block="functions"), backend="cmdstanr", file="fo_result2")

# The software complains about syntax errors

# # the next lines of code show the coupling of brms to the old integrators, works without a problem
# fo_model <- "
#   real[] ode_fo(real t, //time
#   real [] y,         // the rates
#   real [] theta,     // the parameters
#   real [] x_r,       // data constant (not used)
#   int[] x_i){        // data constant (not used)
#   real dydt[1];      // dimension of ODEs
#   dydt[1] = - theta[1] * y[1];                   // first ODE
#    return dydt;        // returns a 3-dimensional array     
#   }
# 
# // this is the function call from brms, integration of ODEs:
# real fo_ode(real t, real A0, real k1) {
#   real y0[1]; //one initial value
#   real theta[1]; //one parameter next to A0
#   real y[1,1]; //ODE solution
#   y0[1]=A0; //initial values
#   theta[1]=k1;
# 
#   y=integrate_ode_rk45(ode_fo, y0,0,rep_array(t,1),theta,rep_array(0.0,0), rep_array(1,1), 0.00001,0.00001,100);
# // Return relevant values
#     return(y[1,1]);
# }
# 
# "
# # t=0 not allowed, will be estimated
# 
# df_sim_ode <- df_sim[-1,]
# 
# fo_formula <- bf(y_obs~fo_ode(time,A0,k1),
#                  A0~1,
#                  k1~1,
#                  nl=TRUE)
# 
# fo_priors <- c(prior(normal(100,10),nlpar=A0),
#                prior(normal(0.2,0.1), nlpar=k1), 
#                prior(cauchy(0,10), class=sigma)
# )
# 
# fo_result1 <- brm(data=df_sim_ode,family = gaussian,formula=fo_formula,prior=fo_priors, inits=0,iter=4000,chains=4, cores=4, stanvars = stanvar(scode=fo_model, block="functions"), file="fo_result1")
