# Q5: grid approximation for ra_prospect (10 grids)
# modified from the class grid approximation example

# rm(list=ls())  # remove all variables

N = 5  # number of subjects
T = 140 # number of trials per subject

# cert[t]: outcome of a certain option on trial t
# gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
# gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.
# 
# evSafe: expected value of a certain (safe) option
# evGamble: expected value of a risky option (gamble)
# pGamble   # probability of choosing a gamble on each trial
# free parameters: rho, tau, lambda

ra_prospect <- function(param, gain, loss, cert, gamble){
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
  evSafe   = cert[t]^rho
  evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
  pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
  pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999

  tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
  sum_minusLL = sum_minusLL + tmp_minusLL
}

  sum_minusLL
}

ra_example <- read.table("ra_exampleData.txt", header = TRUE)

subj_list <- c(2, 3, 4, 6, 7)

grid_mean_10 <- data.frame(
  subjID = subj_list,
  rho = NA,
  lambda = NA,
  tau = NA
)

for(n in 1:N){
  subj_data <- subset(ra_example, subjID==subj_list[n])

# define grid
rho_grid <- seq(0.01, 2, length.out = 10)
lambda_grid <- seq(0.01, 10, length.out = 10)
tau_grid <- seq(0.01, 5, length.out = 10)

grid_results <- expand.grid(
  rho = rho_grid,
  lambda = lambda_grid,
  tau = tau_grid
)

# define prior
prior <- rep(1, nrow(grid_results))

# compute likelihood
grid_results$minusLL <- NA

for (i in 1:nrow(grid_results)) {
  grid_results$minusLL[i] <- ra_prospect(
    param = c(grid_results$rho[i],
              grid_results$lambda[i],
              grid_results$tau[i]),
    gain = subj_data$gain,
    loss = subj_data$loss,
    cert = subj_data$cert,
    gamble = subj_data$gamble
  )
}

# compute likelihood
log_likelihood <- -grid_results$minusLL
likelihood <- exp(log_likelihood - max(log_likelihood))

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior, so it sums to 1
grid_results$posterior <- unstd.posterior / sum(unstd.posterior)

# posterior means
grid_mean_10$rho[n] <- sum(grid_results$rho * grid_results$posterior)
grid_mean_10$lambda[n] <- sum(grid_results$lambda * grid_results$posterior)
grid_mean_10$tau[n] <- sum(grid_results$tau * grid_results$posterior)
}
print(grid_mean_10)

# marginal posterior distributions
rho_post <- aggregate(posterior ~ rho, data = grid_results, sum)
lambda_post <- aggregate(posterior ~ lambda, data = grid_results, sum)
tau_post <- aggregate(posterior ~ tau, data = grid_results, sum)

plot(rho_post$rho, rho_post$posterior, type = "h",
     xlab = "rho", ylab = "Posterior probability")

plot(lambda_post$lambda, lambda_post$posterior, type = "h",
     xlab = "lambda",
     ylab = "Posterior probability")

plot(tau_post$tau, tau_post$posterior, type = "h",
     xlab = "tau",
     ylab = "Posterior probability")
