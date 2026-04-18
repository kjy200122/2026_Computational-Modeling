N = 5  # number of subjects
T = 140 # number of trials per subject

ra_prospect <- function(param, ra_data)  {
  # parameters
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  
  # columns
  cert <- ra_data$cert
  gain <- ra_data$gain
  loss <- ra_data$loss
  gamble <- ra_data$gamble
  
  sum_minusLL <- 0
  
  for (t in 1:T){
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

ra_noLA <- function(param, ra_data)  {
  # parameters
  rho <- param[1]
  lambda <- 1
  tau <- param[2]
  
  # columns
  cert <- ra_data$cert
  gain <- ra_data$gain
  loss <- ra_data$loss
  gamble <- ra_data$gamble
  
  T <- nrow(ra_data)
  sum_minusLL <- 0
  
  for (t in 1:T){
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

ra_noRA <- function(param, ra_data)  {
  # parameters
  rho <- 1
  lambda <- param[1]
  tau <- param[2]
  
  # columns
  cert <- ra_data$cert
  gain <- ra_data$gain
  loss <- ra_data$loss
  gamble <- ra_data$gamble
  
  T <- nrow(ra_data)
  sum_minusLL <- 0
  
  for (t in 1:T){
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

