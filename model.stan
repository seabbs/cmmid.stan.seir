functions {
  real[] seir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      int N = x_i[1];
      real tswitch = x_r[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      real a = theta[3];
      real i0 = theta[4]; // initial number of infected people
      real e0 = theta[5]; // initial number of infected people
      real init[4] = {N - i0 - e0, e0, i0, 0}; // initial values

      real S = y[1] + init[1];
      real E = y[2] + init[2];
      real I = y[3] + init[3];
      real R = y[4] + init[4];
      
      real dS_dt = -beta * I * S / N;
      real dE_dt =  beta * I * S / N - a * E;
      real dI_dt = a * E - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  real t0;
  real tswitch; // date of introduction of control measures
  real ts[n_days];
  int N; // population size
  int cases[n_days];
}
transformed data {
  int x_i[1] = { N }; //formatting to feed the ODE function
  real x_r[1] = {tswitch};
}
parameters {
  real<lower=0> gamma; // SEIR parameters
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> phi_inv; 
  real<lower=0> i0; // number of infected people inititally
  real<lower=0> e0; // number of exposed people inititally
}
transformed parameters{
  real y[n_days, 4];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;
  real theta[5];
  theta = {beta, gamma, a, i0, e0};
  y = integrate_ode_rk45(seir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
  for (i in 1:n_days-1){
    incidence[i] = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]); //-(E(t+1) - E(t) + S(t+1) - S(t))
  }
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  a ~ normal(0.4, 0.5);
  phi_inv ~ exponential(1);
  i0 ~ normal(0, 5);
  e0 ~ normal(0, 5);
  
  // observation model
  cases[1:(n_days-1)] ~ neg_binomial_2(incidence, phi);
}
generated quantities {
  real R0 = beta / gamma;
  real Reff[n_days]; // R0 but taking into account environmental changes
  real recovery_time = 1 / gamma;
  real incubation_time = 1 / a;
  real pred_cases[n_days-1];
  pred_cases = neg_binomial_2_rng(incidence, phi);
  for (i in 1:n_days)
    Reff[i] = beta / gamma;
}