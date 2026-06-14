// Non-hierarchical version: each subject is estimated independently
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
  // Subject-level parameters (no hyperparameters!)
  vector<lower=0,upper=1>[N] alpha;
  vector<lower=0,upper=5>[N] beta;
}

model {
  // Weak priors on individual parameters (no group-level prior)
  alpha ~ uniform(0, 1);
  beta  ~ uniform(0, 5);

  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev;
    real PE;

    ev = initV;

    for (t in 1:(Tsubj[i])) {
      choice[i,t] ~ categorical_logit( beta[i]*ev );

      PE = outcome[i,t] - ev[choice[i,t]];
      ev[choice[i,t]] = ev[choice[i,t]] + alpha[i] * PE;
    }
  }
}
