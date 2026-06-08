data {
  int<lower=1> N; // number of subjects
  int<lower=1> T;
  array[N, T] int<lower=0, upper=1> gamble;
  array[N, T] real cert;
  array[N, T] real<lower=0> gain;
  array[N, T] real<lower=0> loss;  // absolute loss amount
}

parameters {
  vector<lower=0, upper=2>[N] rho;
  vector<lower=0, upper=5>[N] lambda;
  vector<lower=0, upper=10>[N] tau;
}

model {
  // ra_prospect: Original model in Sokol-Hessner et al. 2009 PNAS

  rho    ~ uniform(0, 2);
  lambda ~ uniform(0, 5);
  tau    ~ uniform(0, 10);
  
  for (subj in 1:N) {
    for (t in 1:T) {
      real evSafe;    
      real evGamble;
      real pGamble;
  
      evSafe   = pow(cert[subj, t], rho[subj]);
      evGamble = 0.5 * (pow(gain[subj, t], rho[subj]) -
                  lambda[subj] * pow(loss[subj, t], rho[subj]));
      pGamble  = inv_logit(tau[subj] * (evGamble - evSafe));

      gamble[subj, t] ~ bernoulli(pGamble);
    }
  }
}

generated quantities {
  // For posterior predictive check
  array[N, T] real y_pred;

  for (subj in 1:N) {
    for (t in 1:T) {
      real evSafe;
      real evGamble;
      real pGamble;
  
      evSafe   = pow(cert[subj, t], rho[subj]);
      evGamble = 0.5 * (pow(gain[subj, t], rho[subj]) -
                  lambda[subj] * pow(loss[subj, t], rho[subj]));
      pGamble  = inv_logit(tau[subj] * (evGamble - evSafe));
  
      y_pred[subj, t] = bernoulli_rng(pGamble);
    }
  }
}