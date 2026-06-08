data {
  int<lower=1> N; // number of people
  int<lower=1> T;
  int<lower=0, upper=1> gamble[N, T]; // vector -> matrix
  real cert[N, T];
  real<lower=0> gain[N, T];
  real<lower=0> loss[N, T];  // absolute loss amount
}
transformed data {
}
parameters {
  real<lower=0, upper=2>  rho[N]; // scalar -> vector
  real<lower=0, upper=5> lambda[N];
  real<lower=0, upper=10> tau[N]; 
}
transformed parameters {
}
model {
  // ra_prospect: Original model in Sokol-Hessner et al 2009 PNAS
  // for a single subject
  for (i in 1:N) {  // for loop
    rho[i]    ~ uniform(0, 2);
    lambda[i] ~ uniform(0, 5);
    tau[i]    ~ uniform(0, 10);

    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
      evSafe   = pow(cert[i,t],   rho[i]);     
      evGamble = 0.5 * (pow(gain[i,t], rho[i]) - lambda[i] * pow(loss[i,t], rho[i])); 
      pGamble  = inv_logit(tau[i] * (evGamble - evSafe));  
      gamble[i,t] ~ bernoulli(pGamble);  
    }
  }
}

generated quantities {
  real y_pred[N, T];

  for (i in 1:N) { // for loop
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
      evSafe   = pow(cert[i,t],   rho[i]);
      evGamble = 0.5 * (pow(gain[i,t], rho[i]) - lambda[i] * pow(loss[i,t], rho[i]));
      pGamble  = inv_logit(tau[i] * (evGamble - evSafe));
      y_pred[i,t] = bernoulli_rng(pGamble);
    }
  }
}