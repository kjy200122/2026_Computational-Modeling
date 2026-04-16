N = 5  # number of subjects
T = 140 # number of trials per subject

# for a single subject
mle_ra_prospect <- function(param, T = trial, cert, gain, loss, gamble)  {
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  rho         = param[1]
  lambda      = param[2]
  tau         = param[3]
  
  for (t in 1:T) {
    # cert[t]: outcome of a certain option on trial t
    # gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
    # gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.
    # 
    # evSafe: expected value of a certain (safe) option
    # evGamble: expected value of a risky option (gamble)
    # pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
    
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  sum_minusLL
}

mle_ra_noLA <- function(param, T = trial, cert, gain, loss, gamble)  {
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  rho         = param[1]
  tau         = param[2]
  
  for (t in 1:T) {
    
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  sum_minusLL
}

mle_ra_noRA <- function(param, T = trial, cert, gain, loss, gamble)  {
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  lambda      = param[1]
  tau         = param[2]
  
  for (t in 1:T) {
    
    evSafe   = cert[t]
    evGamble = 0.5*(gain[t]- lambda*abs(loss[t])) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  sum_minusLL
}
