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
# read ra_exampleData
ra_example <- read.table("ra_exampleData.txt", header=TRUE)

subj_list <- c(2, 3, 4, 6, 7)

param1_init <- runif(1)
param2_init <- runif(2)
param3_init <- runif(3)

param_ra_prospect_low <- c(0, 0, 0); param_ra_prospect_up <- c(2, 10, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)

##########################
## MLE                  ##
##########################

# for a single subject (subjectID = 2)
sub2 <- subset(ra_example, subjID==2)

# Call general purpose optimization rountine
mle_ra_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=sub2$gain, loss=sub2$loss, cert=sub2$cert, gamble=sub2$gamble)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_ra_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=sub2$gain, loss=sub2$loss, cert=sub2$cert, gamble=sub2$gamble)

  # Replace the results if the latest optimization yields better result
  if(temp_ra_prospect$value < mle_ra_prospect$value) mle_ra_prospect <- temp_ra_prospect
  }

# Save the MLE parameter estimates
parm_ra_prospect <- mle_ra_prospect$par
value_ra_prospect <- mle_ra_prospect$value

# for all subjects
parm_ra_prospect_all <- data.frame(subjID=subj_list, rho=NA, lambda=NA, tau=NA, minusLL=NA)
for(n in 1:N){
subj_data <- subset(ra_example, subjID==subj_list[n])
param3_init <- runif(3)

mle_ra_prospect_all <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
# Try many different inits to escape from the local maxima
for (i in 1:100)
{
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_ra_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
  
  # Replace the results if the latest optimization yields better result
if(temp_ra_prospect$value < mle_ra_prospect_all$value) mle_ra_prospect_all <- temp_ra_prospect}
# Save the MLE parameter estimates
parm_ra_prospect_all$rho[n] <- mle_ra_prospect_all$par[1]
parm_ra_prospect_all$lambda[n] <- mle_ra_prospect_all$par[2]
parm_ra_prospect_all$tau[n] <- mle_ra_prospect_all$par[3]
parm_ra_prospect_all$minusLL[n] <- mle_ra_prospect_all$value
}
print(round(parm_ra_prospect, 4)) #subject 2 parameters
print(round(value_ra_prospect, 4)) #subject 2 minusLL
print(round(parm_ra_prospect_all, 4)) #all subjects parameters and minusLL

ra_noLA <- function(param, gain, loss, cert, gamble){
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
  
  sum_minusLL
}
param_ra_noLA_low <- c(0, 0); param_ra_noLA_up <- c(2, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)

ra_noRA <- function(param, gain, loss, cert, gamble){
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
  
  sum_minusLL
}
param_ra_noRA_low <- c(0, 0); param_ra_noRA_up <- c(10, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)

# Do the MLE
parm_ra_noLA <- data.frame(subjID=subj_list, rho=NA, lambda=1, tau=NA, minusLL=NA)
for(n in 1:N){
  subj_data <- subset(ra_example, subjID==subj_list[n])
  param2_init <- runif(2)
  
  mle_ra_noLA <- optim(param2_init, ra_noLA, method="L-BFGS-B", lower=param_ra_noLA_low, upper=param_ra_noLA_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
  # Try many different inits to escape from the local maxima
  for (i in 1:100)
  {
    # Re-generate random inits. Is it the best way to do this?
    param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
    
    # Do the MLE again
    temp_ra_noLA <- optim(param2_init, ra_noLA, method="L-BFGS-B", lower=param_ra_noLA_low, upper=param_ra_noLA_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_ra_noLA$value < mle_ra_noLA$value) mle_ra_noLA <- temp_ra_noLA}
  # Save the MLE parameter estimates
  parm_ra_noLA$rho[n] <- mle_ra_noLA$par[1]
  parm_ra_noLA$lambda[n] <- 1
  parm_ra_noLA$tau[n] <- mle_ra_noLA$par[2]
  parm_ra_noLA$minusLL[n] <- mle_ra_noLA$value
}
print(round(parm_ra_noLA, 4)) #all subjects parameters and minusLL

parm_ra_noRA <- data.frame(subjID=subj_list, rho=1, lambda=NA, tau=NA, minusLL=NA)
for(n in 1:N){
  subj_data <- subset(ra_example, subjID==subj_list[n])
  param2_init <- runif(2)
  
  mle_ra_noRA <- optim(param2_init, ra_noRA, method="L-BFGS-B", lower=param_ra_noRA_low, upper=param_ra_noRA_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
  # Try many different inits to escape from the local maxima
  for (i in 1:100)
  {
    # Re-generate random inits. Is it the best way to do this?
    param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
    
    # Do the MLE again
    temp_ra_noRA <- optim(param2_init, ra_noRA, method="L-BFGS-B", lower=param_ra_noRA_low, upper=param_ra_noRA_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_ra_noRA$value < mle_ra_noRA$value) mle_ra_noRA <- temp_ra_noRA}
  # Save the MLE parameter estimates
  parm_ra_noRA$rho[n] <- 1
  parm_ra_noRA$lambda[n] <- mle_ra_noRA$par[1]
  parm_ra_noRA$tau[n] <- mle_ra_noRA$par[2]
  parm_ra_noRA$minusLL[n] <- mle_ra_noRA$value
}
print(round(parm_ra_noRA, 4)) #all subjects parameters and minusLL

# Calculate AIC
k_ra_prospect <- 3
k_ra_noLA <- 2
k_ra_noRA <- 2

AIC_ra_prospect <- 2*k_ra_prospect + 2*parm_ra_prospect_all$minusLL
AIC_ra_noLA <- 2*k_ra_noLA + 2*parm_ra_noLA$minusLL
AIC_ra_noRA <- 2*k_ra_noRA + 2*parm_ra_noRA$minusLL
sum_AIC_ra_prospect <- sum(AIC_ra_prospect)
sum_AIC_ra_noLA <- sum(AIC_ra_noLA)
sum_AIC_ra_noRA <- sum(AIC_ra_noRA)

# Calculate BIC
BIC_ra_prospect <- log(T)*k_ra_prospect + 2*parm_ra_prospect_all$minusLL
BIC_ra_noLA <- log(T)*k_ra_noLA + 2*parm_ra_noLA$minusLL
BIC_ra_noRA <- log(T)*k_ra_noRA + 2*parm_ra_noRA$minusLL
sum_BIC_ra_prospect <- sum(BIC_ra_prospect)
sum_BIC_ra_noLA <- sum(BIC_ra_noLA)
sum_BIC_ra_noRA <- sum(BIC_ra_noRA)

#Model Comparison
model_comp <- data.frame(
  Model = c("ra_prospect", "ra_noLA", "ra_noRA"),
  AIC = c(sum_AIC_ra_prospect, sum_AIC_ra_noLA, sum_AIC_ra_noRA),
  BIC = c(sum_BIC_ra_prospect, sum_BIC_ra_noLA, sum_BIC_ra_noRA)
)

print(model_comp)
