data {
  int<lower=1> N;                          // number of subjects
  int<lower=1> T;                          // number of trials per subject
  int<lower=0, upper=1> gamble[N, T];
  real cert[N, T];
  real<lower=0> gain[N, T];
  real<lower=0> loss[N, T];                // absolute loss amount
}
transformed data {
}
parameters {
  real<lower=0, upper=2>  rho[N];
  real<lower=0, upper=5>  lambda[N];
  real<lower=0, upper=10> tau[N];
}
transformed parameters {
}
model {
  // ra_prospect: Original model in Sokol-Hessner et al 2009 PNAS
  // modified to fit all N subjects with a for loop within Stan
  for (i in 1:N) {

    rho[i]    ~ uniform(0, 2);
    lambda[i] ~ uniform(0, 5);
    tau[i]    ~ uniform(0, 10);

    for (t in 1:T) {

      real evSafe;
      real evGamble;
      real pGamble;
      // loss[i,t] = absolute amount of loss (pre-converted in R)
      evSafe   = pow(cert[i, t], rho[i]);
      evGamble = 0.5 * (pow(gain[i, t], rho[i]) - lambda[i] * pow(loss[i, t], rho[i]));
      pGamble  = inv_logit(tau[i] * (evGamble - evSafe));
      gamble[i, t] ~ bernoulli(pGamble);
    }
  }
}
generated quantities {
  // For posterior predictive check
  real y_pred[N, T];
  for (i in 1:N) {
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
      evSafe   = pow(cert[i, t], rho[i]);
      evGamble = 0.5 * (pow(gain[i, t], rho[i]) - lambda[i] * pow(loss[i, t], rho[i]));
      pGamble  = inv_logit(tau[i] * (evGamble - evSafe));
      // generate posterior prediction for current trial
      y_pred[i, t] = bernoulli_rng(pGamble);
    }
  }
}
