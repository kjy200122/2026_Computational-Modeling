// Non-hierarchical model with alpha and beta
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1,upper=T> Tsubj[N];
  int<lower=-1,upper=2> choice[N,T];
  real outcome[N,T];
}

transformed data {
  vector[2] initV;
  initV = rep_vector(0.0, 2);
}

parameters {
  // Subject-level parameters - estimated independently, no group structure
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0, upper=5>[N] beta;
}
model {
  // Weakly informative priors directly on each subject's parameters
  alpha ~ uniform(0, 1);   // or beta(1, 1), same thing
  beta  ~ uniform(0, 5);   // or normal(0, 5) truncated, etc.

  // Subject loop and trial loop - identical to before
  for (i in 1:N) {
    vector[2] ev;
    real PE;
    ev = initV;
    for (t in 1:(Tsubj[i])) {
      choice[i,t] ~ categorical_logit( beta[i] * ev );
      PE = outcome[i,t] - ev[choice[i,t]];
      ev[choice[i,t]] = ev[choice[i,t]] + alpha[i] * PE;
    }
  }
}
