data {
  int<lower=1> N; // to estimate all subjects
  int<lower=1> T;
  int<lower=0, upper=1> gamble[N, T]; 
  matrix[N, T] cert; // use matrix to apply to all subjects
  matrix<lower=0>[N, T] gain;
  matrix<lower=0>[N, T] loss;  // absolute loss amount
}
transformed data {
}
parameters {
  vector<lower=0, upper=2>[N] rho; // set parameters as vector
  vector<lower=0, upper=5>[N] lambda;
  vector<lower=0, upper=10>[N] tau; 
 }

transformed parameters {
}
model {
  // ra_prospect: Original model in Sokol-Hessner et al 2009 PNAS
  // for all subjects
  rho    ~ uniform(0, 2);
  lambda ~ uniform(0, 5);
  tau    ~ uniform(0, 10);
  
for (n in 1:N) { // use for loop to estimate all subjects
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
  int y_pred[N, T];

for (n in 1:N) {
  for (t in 1:T) {
    real evSafe;
    real evGamble;
    real pGamble;

    evSafe     = pow(cert[n, t], rho[n]);
    evGamble   = 0.5 * (pow(gain[n, t], rho[n]) - lambda[n] * pow(loss[n, t], rho[n]));
    pGamble    = inv_logit(tau[n] * (evGamble - evSafe));

    // generate posterior prediction for current trial
    y_pred[n, t] = bernoulli_rng(pGamble);
  }
 }
}
