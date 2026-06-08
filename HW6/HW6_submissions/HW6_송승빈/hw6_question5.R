N = 5  # number of subjects
T = 140 # number of trials per subject

# dataset
data <- read.table("ra_exampleData.txt", header = TRUE)

# model (ra_prospect)
ra_prospect <- function(param, cert, gain, loss, gamble, T){
  loglikelihood = 0  
  
  for (t in 1:T) {
    # cert[t]: outcome of a certain option on trial t
    # gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
    # gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.
    # 
    # evSafe: expected value of a certain (safe) option
    # evGamble: expected value of a risky option (gamble)
    # pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
    
    evSafe   = cert[t]^param[1]
    evGamble = 0.5*(gain[t]^param[1] - param[2]*abs(loss[t])^param[1]) 
    pGamble  = 1 / (1 + exp(param[3]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_loglik = log(pGamble)*gamble[t] + log(1-pGamble)*(1-gamble[t])  # likelihood of trial t (instead of minus log likelihood)
    loglikelihood = loglikelihood + tmp_loglik
  }
  
  loglikelihood
}

# define grid
rho_grid <- seq( from=0.01, to=2, length.out=30 )
lambda_grid <- seq( from=0, to=5, length.out=30 )
tau_grid <- seq( from=0, to=10, length.out=30 )


# prior (uniform prior)
prior <- array(1, dim=c(30,30,30))

posterior <- vector("list", N)

for (s in 1:N){
  
  curr_data <- data[data$subjID == s, ]
  
  loglikelihood <- array(NA, dim=c(30,30,30))

  for (i in 1:10){
    for (j in 1:10){
      for (k in 1:10){
        param <- c(
          rho_grid[i],
          lambda_grid[j],
          tau_grid[k]
        )
        # likelihood
        loglikelihood[i,j,k] <-
          ra_prospect(
              param,
              curr_data$cert,
              curr_data$gain,
              curr_data$loss,
              curr_data$gamble,
              T
          )
      }
    }
  }

  # product of likelihood and prior(which is just 1)
  loglikelihood <- loglikelihood - max(loglikelihood, na.rm = TRUE) # handling error
  loglikelihood[is.na(loglikelihood)] <- -1e10 # handling error
  
  unstd.posterior <- exp(loglikelihood + log(1))
  unstd.posterior <- unstd.posterior + 1e-10

  # standardize the posterior, so it sums to 1
  posterior[[s]] <- unstd.posterior / sum(unstd.posterior)
}

summary(data$gain)
summary(data$loss)
summary(data$cert)

# ========= Plot results ======================================================

# for subject 1
post <- posterior[[1]]

# rho에 대한 marginal posterior
rho_post <- apply(post, 1, sum)  
summary(rho_post)

summary(loglikelihood)
range(loglikelihood)

# lambda
lambda_post <- apply(post, 2, sum)  
summary(lambda_post)

# tau
tau_post <- apply(post, 3, sum)  
summary(tau_post)

# plot
plot(rho_grid, rho_post, type="l", lwd=2, xlab="rho", ylab="Posterior", main="Posterior of rho")
