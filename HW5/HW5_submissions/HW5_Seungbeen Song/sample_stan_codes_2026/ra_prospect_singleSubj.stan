data {
  int<lower=1> N;
  int<lower=1> T;
  
  int<lower=0, upper=1> gamble[N,T]; 
  matrix<lower=0>[N,T] cert;
  matrix<lower=0>[N,T] gain;
  matrix<lower=0>[N,T] loss;  // absolute loss amount
}
transformed data {
}
parameters {
  vector<lower=0, upper=2>[N] rho;
  vector<lower=0, upper=5>[N] lambda;
  vector<lower=0, upper=10>[N] tau; 
}
transformed parameters {
}
model {
  // ra_prospect: Original model in Sokol-Hessner et al 2009 PNAS
  for (s in 1:N){
  rho[s]    ~ uniform(0, 2);
  lambda[s] ~ uniform(0, 5);
  tau[s]   ~ uniform(0, 10);
  
  for (t in 1:T) {
    
    real evSafe;    
    real evGamble;
    real pGamble;

    // loss[t]=absolute amount of loss (pre-converted in R)
    evSafe   = pow(cert[s,t], rho[s]);
    evGamble = 0.5 * (pow(gain[s,t], rho[s]) - lambda[s] * pow(loss[s,t], rho[s]));
    pGamble  = inv_logit(tau[s] * (evGamble - evSafe));
    gamble[s,t] ~ bernoulli(pGamble);
  }
 }
}
generated quantities{
  // For posterior predictive check
  real y_pred[N,T];
  
  for (s in 1:N){
  for (t in 1:T) {
    real evSafe;
    real evGamble;
    real pGamble;

    evSafe     = pow(cert[s,t], rho[s]);
    evGamble   = 0.5 * (pow(gain[s,t], rho[s]) - lambda[s] * pow(loss[s,t], rho[s]));
    pGamble    = inv_logit(tau[s] * (evGamble - evSafe));

    // generate posterior prediction for current trial
    y_pred[s,t] = bernoulli_rng(pGamble);
  }
 }
}
