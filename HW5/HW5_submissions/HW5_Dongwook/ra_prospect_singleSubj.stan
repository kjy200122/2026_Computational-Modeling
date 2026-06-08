data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=0, upper=1> gamble[N, T]; 
  real cert[N, T];
  real<lower=0> gain[N, T];
  real<lower=0> loss[N, T];  // absolute loss amount
}

parameters {
  vector<lower=0, upper=2>[N] rho;
  vector<lower=0, upper=5>[N] lambda;
  vector<lower=0, upper=10>[N] tau; 
}

model {
  // ra_prospect: Original model in Sokol-Hessner et al 2009 PNAS
  // for a single subject
  rho    ~ uniform(0, 2);
  lambda ~ uniform(0, 5);
  tau    ~ uniform(0, 10);
  for (n in 1:N) {
    for (t in 1:T) {
      
      real evSafe;    
      real evGamble;
      real pGamble;
  
      // loss[t]=absolute amount of loss (pre-converted in R)
      evSafe   = pow(cert[n, t], rho[n]);
      evGamble = 0.5 * (pow(gain[n, t], rho[n]) - lambda[n] * pow(loss[n, t], rho[n]));
      pGamble  = inv_logit(tau[n] * (evGamble - evSafe));
      gamble[n, t] ~ bernoulli(pGamble);
    }
  }
}
generated quantities{
  // For posterior predictive check
  real y_pred[T];
  for (n in 1:N) {
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
  
      evSafe     = pow(cert[n, t], rho[n]);
      evGamble   = 0.5 * (pow(gain[n, t], rho[n]) - lambda[n] * pow(loss[n, t], rho[n]));
      pGamble    = inv_logit(tau[n] * (evGamble - evSafe));
  
      // generate posterior prediction for current trial
      y_pred[t] = bernoulli_rng(pGamble);
    }
  }
}
