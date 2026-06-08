data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=0, upper=1> gamble[N, T];
  real cert[N, T];
  real<lower=0> gain[N, T];
  real<lower=0> loss[N, T];  // absolute loss amount
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
  // for a single subject
  rho ~ uniform(0, 2);
  lambda ~ uniform(0, 5);
  tau ~ uniform(0, 10);

  for (s in 1:N) {
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
      
      // loss[t]=absolute amount of loss (pre-converted in R)
      evSafe   = pow(cert[s, t], rho[s]);
      evGamble = 0.5 * (pow(gain[s, t], rho[s]) - lambda[s] * pow(loss[s, t], rho[s]));
      pGamble  = inv_logit(tau[s] * (evGamble - evSafe));
      gamble[s, t] ~ bernoulli(pGamble);
    }
  }
}
generated quantities{
  // For posterior predictive check
  array[N, T] real y_pred;

  for (s in 1:N) {
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;

      evSafe   = pow(cert[s, t], rho[s]);
      evGamble = 0.5 * (pow(gain[s, t], rho[s]) - lambda[s] * pow(loss[s, t], rho[s]));
      pGamble  = inv_logit(tau[s] * (evGamble - evSafe));

      // generate posterior prediction for current trial
      y_pred[s, t] = bernoulli_rng(pGamble);
    }
  }
}
