N = 5  # number of subjects
T = 140 # number of trials per subject

# ================== Problem 1 & 2 ================== 

param_init <-runif(3) # initial values
param_low <- c(0,0,0); param_up <- c(2,10,5) # parameter bounds

### For a single subject

data <- read.table("ra_exampleData.txt", header = TRUE)
View(data)

data_sub2 <- data[data$subjID == 2, ]

ra_prospect <- function (param, cert, gain, loss, gamble, T){
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
  
    evSafe   = cert[t]^param[1]
    evGamble = 0.5*(gain[t]^param[1] - param[2]*abs(loss[t])^param[1]) 
    pGamble  = 1 / (1 + exp(param[3]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  sum_minusLL
}

# Run MLE
ra_pro_mle <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=data_sub2$cert, gain=data_sub2$gain, loss=data_sub2$loss, gamble=data_sub2$gamble, T=140)

# Try many different initial values
for (i in 1:100){
  param_init <- runif(3)
  ra_pro_temp <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=data_sub2$cert, gain=data_sub2$gain, loss=data_sub2$loss, gamble=data_sub2$gamble, T=140)
  
  if (ra_pro_temp$value < ra_pro_mle$value) ra_pro_mle <- ra_pro_temp
}

print(ra_pro_mle$par)


### For all subjects 

all_parameters <- data.frame(
  subjID = unique(data$subjID),
  rho = numeric(N),
  lambda = numeric(N),
  tau = numeric(N)
)

idx = 1

for (i in unique(data$subjID)){
  curr_data <- data[data$subjID == i, ]
  param_init <- runif(3)
  
  # Run MLE
  rp_mle <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
  
  # Try many different initial values
  for (j in 1:100){
    param_init <- runif(3)
    
    rp_temp <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
    
    if (rp_temp$value < rp_mle$value) rp_mle <- rp_temp
  }
  
  all_parameters$rho[idx] <- rp_mle$par[1]
  all_parameters$lambda[idx] <- rp_mle$par[2]
  all_parameters$tau[idx] <- rp_mle$par[3]
  
  idx = idx + 1
  
}

print(all_parameters)


# ================== Problem 3 ================== 

### Model without lambda

ra_noLA <- function (param, cert, gain, loss, gamble, T){
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    # free parameters: rho, tau
    
    evSafe   = cert[t]^param[1]
    evGamble = 0.5*(gain[t]^param[1] - abs(loss[t])^param[1])
    pGamble  = 1 / (1 + exp(param[2]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  sum_minusLL
}

param_low <- c(0,0); param_up <- c(2,5) # parameter bounds

# Create table for all participants with parameter values, AIC, and BIC
df_ra_noLA <- data.frame(
  subjID = unique(data$subjID),
  rho = numeric(N),
  tau = numeric(N),
  AIC = numeric(N),
  BIC = numeric(N)
)

idx = 1

for (i in unique(data$subjID)){
  curr_data <- data[data$subjID == i, ]
  param_init <- runif(2)
  
  # Run MLE
  rnl_mle <- optim(param_init, ra_noLA, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
  
  # Try many different initial values
  for (j in 1:100){
    param_init <- runif(2)
    
    rnl_temp <- optim(param_init, ra_noLA, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
    
    if (rnl_temp$value < rnl_mle$value) rnl_mle <- rnl_temp
  }
  
  k_rnl <- 2
  num <- nrow(curr_data)
  
  idv_AIC_rnl <- 2*rnl_mle$value + 2*k_rnl 
  idv_BIC_rnl <- 2*rnl_mle$value + k_rnl*log(num)

  
  df_ra_noLA$rho[idx] <- rnl_mle$par[1]
  df_ra_noLA$tau[idx] <- rnl_mle$par[2]
  df_ra_noLA$AIC[idx] <- idv_AIC_rnl
  df_ra_noLA$BIC[idx] <- idv_BIC_rnl
  
  idx = idx + 1

}

AIC_noLA = sum(df_ra_noLA$AIC)
BIC_noLA = sum(df_ra_noLA$BIC)

print(AIC_noLA)
print(BIC_noLA)



### Model without rho

ra_noRA <- function (param, cert, gain, loss, gamble, T){
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    # free parameters: tau, lambda
    
    evSafe   = cert[t]
    evGamble = 0.5*(gain[t] - param[1]*abs(loss[t])) 
    pGamble  = 1 / (1 + exp(param[2]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  sum_minusLL
}

param_low <- c(0,0); param_up <- c(10,5) # parameter bounds

# Create table for all participants with parameter values, AIC, and BIC
df_ra_noRA <- data.frame(
  subjID = unique(data$subjID),
  lambda = numeric(N),
  tau = numeric(N),
  AIC = numeric(N),
  BIC = numeric(N)
)

idx = 1

for (i in unique(data$subjID)){
  curr_data <- data[data$subjID == i, ]
  param_init <- runif(2)
  
  # Run MLE
  rnr_mle <- optim(param_init, ra_noRA, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
  
  # Try many different initial values
  for (j in 1:100){
    param_init <- runif(2)
    
    rnr_temp <- optim(param_init, ra_noRA, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
    
    if (rnr_temp$value < rnr_mle$value) rnr_mle <- rnr_temp
  }
  
  k_rnr <- 2
  num <- nrow(curr_data)
  
  idv_AIC_rnr <- 2*rnr_mle$value + 2*k_rnr 
  idv_BIC_rnr <- 2*rnr_mle$value + k_rnr*log(num)
  
  df_ra_noRA$lambda[idx] <- rnr_mle$par[1]
  df_ra_noRA$tau[idx] <- rnr_mle$par[2]
  df_ra_noRA$AIC[idx] <- idv_AIC_rnr
  df_ra_noRA$BIC[idx] <- idv_BIC_rnr
  
  idx = idx + 1
  
}

AIC_noRA = sum(df_ra_noRA$AIC)
BIC_noRA = sum(df_ra_noRA$BIC)

print(AIC_noRA)
print(BIC_noRA)


### computing the AIC/BIC of ra_prospect for Problem 3. 3)

AIC_prospect = 0
BIC_prospect = 0

param_low <- c(0,0,0); param_up <- c(2,10,5)

for (i in unique(data$subjID)){
  curr_data <- data[data$subjID == i, ]
  param_init <- runif(3)
  
  # Run MLE
  rp_mle <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
  
  # Try many different initial values
  for (j in 1:100){
    param_init <- runif(3)
    
    rp_temp <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
    
    if (rp_temp$value < rp_mle$value) rp_mle <- rp_temp
  }
  
  k_rp <- 3
  num <- nrow(curr_data)
  
  idv_AIC_rp <- 2*rp_mle$value + 2*k_rp 
  idv_BIC_rp <- 2*rp_mle$value + k_rp*log(num)
  
  AIC_prospect = AIC_prospect + idv_AIC_rp
  BIC_prospect = BIC_prospect + idv_BIC_rp
  
}

print(AIC_prospect)
print(BIC_prospect)
