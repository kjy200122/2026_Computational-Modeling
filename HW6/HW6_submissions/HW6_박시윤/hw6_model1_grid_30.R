dat <- read.table("ra_exampleData.txt", header = TRUE)

ra_prospect <- function(param, gain, loss, cert, gamble) {
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  sum_LL = 0
  
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001
    tmp_LL = log(pGamble)*gamble[t] + log(1-pGamble)*(1-gamble[t])
    sum_LL = sum_LL + tmp_LL
  }
  
  return(sum_LL)
}

# 30 grids
rho_grid <- seq(0.01, 2, length.out = 30)
lambda_grid <- seq(0.01, 10, length.out = 30)
tau_grid <- seq(0.01, 5, length.out = 30)

sub_summary <- data.frame()

for (s in c(2, 3, 4, 6, 7)) {
  sub <- dat[dat$subjID == s, ]
  grid_summary <- expand.grid(
    rho = rho_grid,
    lambda = lambda_grid,
    tau = tau_grid
  )
  grid_summary$LL <- NA

  for (i in 1:nrow(grid_summary)) {
    param <- c(
      grid_summary$rho[i],
      grid_summary$lambda[i],
      grid_summary$tau[i]
    )

    grid_summary$LL[i] <- ra_prospect(
      param,
      gain = sub$gain,
      loss = sub$loss,
      cert = sub$cert,
      gamble = sub$gamble
    )
  }

  # posterior
  grid_summary$posterior <- exp(grid_summary$LL - max(grid_summary$LL))
  grid_summary$posterior <- grid_summary$posterior / sum(grid_summary$posterior)

  # marginal posterior distributions
  rho_post <- aggregate(posterior ~ rho, grid_summary, sum)
  lambda_post <- aggregate(posterior ~ lambda, grid_summary, sum)
  tau_post <- aggregate(posterior ~ tau, grid_summary, sum)

  # plots
  par(mfrow = c(1,3))

  plot(rho_post$rho,
       rho_post$posterior,
       type = "b",
       main = paste("Subject", s, ": rho"))

  plot(lambda_post$lambda,
       lambda_post$posterior,
       type = "b",
       main = paste("Subject", s, ": lambda"))

  plot(tau_post$tau,
       tau_post$posterior,
       type = "b",
       main = paste("Subject", s, ": tau"))
}
