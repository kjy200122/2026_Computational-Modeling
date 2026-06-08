library(ggplot2)
library(reshape2)


data <- read.table("ra_exampleData.txt", header=TRUE, sep="\t")
subj1 <- subset(data, subjID == unique(data$subjID)[1])

# Grid
rho_grid    <- seq(0.1, 2.0, length.out = 10)
lambda_grid <- seq(0.1, 5.0, length.out = 10)
tau_grid    <- seq(0.1, 5.0, length.out = 10)

# posterior
posterior <- array(0, dim = c(10, 10, 10))

# likelihood
log_likelihood <- function(rho, lambda, tau, ra_data) {
  T   <- nrow(ra_data)
  cert   <- ra_data$cert
  gain   <- ra_data$gain
  loss   <- ra_data$loss
  gamble <- ra_data$gamble
  
  loglik <- 0
  for (t in 1:T) {
    evSafe   <- cert[t]^rho
    evGamble <- 0.5 * (gain[t]^rho - lambda * abs(loss[t])^rho)
    pGamble  <- 1 / (1 + exp(tau * (evSafe - evGamble)))
    pGamble  <- pGamble * 0.9998 + 0.0001
    loglik   <- loglik + gamble[t]*log(pGamble) + (1-gamble[t])*log(1-pGamble)
  }
  return(loglik)
}

# Grid calculation
for (i in 1:10) {
  for (j in 1:10) {
    for (k in 1:10) {
      ll <- log_likelihood(rho_grid[i], lambda_grid[j], tau_grid[k], subj1)
      # uniform prior → log_prior = 0
      posterior[i, j, k] <- exp(ll)
    }
  }
}

# normalization
posterior <- posterior / sum(posterior)

# marginal posterior
post_rho    <- apply(posterior, 1, sum)
post_lambda <- apply(posterior, 2, sum)
post_tau    <- apply(posterior, 3, sum)

# draw plots
par(mfrow = c(1, 3))
plot(rho_grid,    post_rho,    type="b", xlab="rho",    ylab="Posterior", main="Marginal: rho")
plot(lambda_grid, post_lambda, type="b", xlab="lambda", ylab="Posterior", main="Marginal: lambda")
plot(tau_grid,    post_tau,    type="b", xlab="tau",    ylab="Posterior", main="Marginal: tau")
par(mfrow = c(1, 1))