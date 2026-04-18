N = 5  # number of subjects
T = 140 # number of trials per subject

# for a single subject

dat <- read.table("ra_exampleData.txt", header = TRUE)

ra_prospect <- function(param, gain, loss, cert, gamble) {
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
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
  
  return(sum_minusLL)
}

ra_noLA <- function(param, gain, loss, cert, gamble) {
  rho <- param[1]
  lambda <- 1
  tau <- param[2]
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  return(sum_minusLL)
}

ra_noRA <- function(param, gain, loss, cert, gamble) {
  rho <- 1
  lambda <- param[1]
  tau <- param[2]
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  return(sum_minusLL)
}

sub2 <- dat[dat$subjID == 2,]

fit_sub2 <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                  lower = c(0, 0, 0), upper = c(2, 10, 5),
                  gain = sub2$gain, loss = sub2$loss,
                  cert = sub2$cert, gamble = sub2$gamble)

for (i in 1:100) {
  temp_fit_sub2 <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                         lower = c(0, 0, 0), upper = c(2, 10, 5),
                         gain = sub2$gain, loss = sub2$loss,
                         cert = sub2$cert, gamble = sub2$gamble)
  
  if (temp_fit_sub2$value < fit_sub2$value) fit_sub2 <- temp_fit_sub2
}

sub2_summary <- data.frame(rho = fit_sub2$par[1],
                           lambda = fit_sub2$par[2],
                           tau = fit_sub2$par[3])

sub_summary <- data.frame()

for (s in c(2, 3, 4, 6, 7)) {
  sub <- dat[dat$subjID == s, ]
  fit_sub <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                   lower = c(0, 0, 0), upper = c(2, 10, 5),
                   gain = sub$gain, loss = sub$loss,
                   cert = sub$cert, gamble = sub$gamble)
  
  for (i in 1:100) {
    temp_fit_sub <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                          lower = c(0, 0, 0), upper = c(2, 10, 5),
                          gain = sub$gain, loss = sub$loss,
                          cert = sub$cert, gamble = sub$gamble)
    
    if (temp_fit_sub$value < fit_sub$value) fit_sub <- temp_fit_sub
  }
  
  sub_summary <- rbind(sub_summary, data.frame(
    subjID = s,
    rho = fit_sub$par[1],
    lambda = fit_sub$par[2],
    tau = fit_sub$par[3])
  )
}

noLA_summary <- data.frame()
noRA_summary <- data.frame()

sum_minusLL_prospect <- 0
sum_minusLL_noLA <- 0
sum_minusLL_noRA <- 0

for (s in c(2, 3, 4, 6, 7)) {
  sub <- dat[dat$subjID == s, ]
  
  fit_pros <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                    lower = c(0, 0, 0), upper = c(2, 10, 5),
                    gain = sub$gain, loss = sub$loss,
                    cert = sub$cert, gamble = sub$gamble)
  
  fit_noLA <- optim(runif(2), ra_noLA, method = "L-BFGS-B",
                    lower = c(0, 0), upper = c(2, 5),
                    gain = sub$gain, loss = sub$loss,
                    cert = sub$cert, gamble = sub$gamble)
  
  fit_noRA <- optim(runif(2), ra_noRA, method = "L-BFGS-B",
                    lower = c(0, 0), upper = c(10, 5),
                    gain = sub$gain, loss = sub$loss,
                    cert = sub$cert, gamble = sub$gamble)
  
  for (i in 1:100) {
    temp_fit_pros <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                           lower = c(0, 0, 0), upper = c(2, 10, 5),
                           gain = sub$gain, loss = sub$loss,
                           cert = sub$cert, gamble = sub$gamble)
  
    temp_fit_noLA <- optim(runif(2), ra_noLA, method = "L-BFGS-B",
                           lower = c(0, 0), upper = c(2, 5),
                           gain = sub$gain, loss = sub$loss,
                           cert = sub$cert, gamble = sub$gamble)
    
    temp_fit_noRA <- optim(runif(2), ra_noRA, method = "L-BFGS-B",
                           lower = c(0, 0), upper = c(10, 5),
                           gain = sub$gain, loss = sub$loss,
                           cert = sub$cert, gamble = sub$gamble)
    
    if (temp_fit_pros$value < fit_pros$value) fit_pros <- temp_fit_pros
    if (temp_fit_noLA$value < fit_noLA$value) fit_noLA <- temp_fit_noLA
    if (temp_fit_noRA$value < fit_noRA$value) fit_noRA <- temp_fit_noRA
  }
  
  sum_minusLL_prospect <- sum_minusLL_prospect + fit_pros$value
  sum_minusLL_noLA <- sum_minusLL_noLA + fit_noLA$value
  sum_minusLL_noRA <- sum_minusLL_noRA + fit_noRA$value
}  

k_prospect <- 3
k_noLA <- 2
k_noRA <- 2

AIC_prospect <- 2*sum_minusLL_prospect + 2*k_prospect*N
AIC_noLA <- 2*sum_minusLL_noLA + 2*k_noLA*N
AIC_noRA <- 2*sum_minusLL_noRA + 2*k_noRA*N

BIC_prospect <- 2*sum_minusLL_prospect + k_prospect*log(T)*N
BIC_noLA <- 2*sum_minusLL_noLA + k_noLA*log(T)*N
BIC_noRA <- 2*sum_minusLL_noRA + k_noRA*log(T)*N

comparison_summary <- data.frame(
  model = c("ra_prospect", "ra_noLA", "ra_noRA"),
  AIC = c(AIC_prospect, AIC_noLA, AIC_noRA),
  BIC = c(BIC_prospect, BIC_noLA, BIC_noRA)
)

print(sub2_summary)
print(sub_summary)
print(comparison_summary)
