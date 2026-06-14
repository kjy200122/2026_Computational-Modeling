// Model with two learning rates (alpha_pos, alpha_neg) and beta
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
  // Hyper(group)-parameters: now 3 parameters
  vector[3] mu_p;
  vector<lower=0>[3] sigma;

  // Subject-level raw parameters
  vector[N] alpha_pos_pr;
  vector[N] alpha_neg_pr;
  vector[N] beta_pr;
}

transformed parameters {
  vector<lower=0,upper=1>[N] alpha_pos;
  vector<lower=0,upper=1>[N] alpha_neg;
  vector<lower=0,upper=5>[N] beta;

  for (i in 1:N) {
    alpha_pos[i] = Phi_approx( mu_p[1] + sigma[1] * alpha_pos_pr[i] );
    alpha_neg[i] = Phi_approx( mu_p[2] + sigma[2] * alpha_neg_pr[i] );
    beta[i]      = Phi_approx( mu_p[3] + sigma[3] * beta_pr[i] ) * 5;
  }
}

model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ normal(0, 1);

  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev;
    real PE;

    ev = initV;

    for (t in 1:(Tsubj[i])) {
      choice[i,t] ~ categorical_logit( beta[i]*ev );

      PE = outcome[i,t] - ev[choice[i,t]];

      // value updating with two learning rates
      if (PE >= 0) {
        ev[choice[i,t]] = ev[choice[i,t]] + alpha_pos[i] * PE;
      } else {
        ev[choice[i,t]] = ev[choice[i,t]] + alpha_neg[i] * PE;
      }
    }
  }
}
