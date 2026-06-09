// model with alpha and beta
data {
  int<lower=1> N;
  int<lower=1> T;               
  int<lower=1,upper=T> Tsubj[N];                 
  int<lower=-1,upper=2> choice[N,T];     
  real outcome[N,T];  // no lower and upper bounds   
}

transformed data {
  vector[2] initV;  // initial values for EV
  initV = rep_vector(0.0, 2);
}

parameters {
  // Subject-level raw parameters (for Matt trick)
  vector<lower=0,upper=1>[N] alpha; // learning rate [0, 1]
  vector<lower=0,upper=5>[N] beta;  // inverse temperature [0, 5]
}

transformed parameters {
}

model {
  // individual parameters
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  
  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev; // expected value
    real PE;      // prediction error

    ev = initV;

    for (t in 1:(Tsubj[i])) {        
      // compute action probabilities
      choice[i,t] ~ categorical_logit( beta[i]*ev );

      // prediction error 
      PE = outcome[i,t] - ev[choice[i,t]];

      // value updating (learning) 
      ev[choice[i,t]] = ev[choice[i,t]] + alpha[i] * PE; 
    }
  }
}
generated quantities {
}
