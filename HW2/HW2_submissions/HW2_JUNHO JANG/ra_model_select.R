N = 5  # number of subjects
T = 140 # number of trials per subject / nrow(ra_data)
source("ra_MLE.R")   # source ra_MLe.R code you programmed for HW2

# data loading
ra_data <- read.table("ra_exampleData.txt", header=TRUE)

#================================== ra_prospect for sub 2  /// HW 2-1========================

# only ID = 2, # for a single subject
data_sub2 <- subset(ra_data, subjID == 2)

#parameters 
param_init <- runif(3)

#bounds
param_prospect_low <- c(0, 0, 0); param_prospect_up <- c(2, 10, 5) #rho, lambda, tau


# Call general purpose optimization rountine
mle_model_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub2)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  
  param_init <- runif(3)
  
  # Do the MLE again
  temp_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub2)
  
  # Replace the results if the latest optimization yields better result
  if(temp_prospect$value < mle_model_prospect$value)  mle_model_prospect <- temp_prospect
  
}

# Save the MLE parameter estimates
parm_prospect <- mle_model_prospect$par
print(parm_prospect)


#==================for all of subjects   /  HW 2-2======================
# checking subjID
subject_id <- unique(ra_data$subjID) #use unique function


# make empty tables
results <- data.frame(subjID = subject_id, rho = NA, lambda = NA, tau = NA)

for (sub in 1:N){
  sub_id <- subject_id[sub]
  
  data_sub <- subset(ra_data, subjID == sub_id)
  
  mle_model_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub)
  
  for (i in 1:100) {
    
    param_init <- runif(3)
    
    # Do the MLE again
    temp_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub)
    
    # Replace the results if the latest optimization yields better result
    if(temp_prospect$value < mle_model_prospect$value)  mle_model_prospect <- temp_prospect
    
  }
  results$rho[sub] <- mle_model_prospect$par[1]
  results$lambda[sub] <- mle_model_prospect$par[2]
  results$tau[sub] <- mle_model_prospect$par[3]
}

print(results)


#================================== ra_noLA/// HW 3-1========================
results_noLA <- data.frame(subjID = subject_id, AIC = NA, BIC = NA)

total_noLA <- 0 # for all subjects AIC, BIC
k_noLA <- 2 #every subjects has two parameters


#parameters 
param_init <- runif(2)

#bounds
param_noLA_low <- c(0, 0); param_noLA_up <- c(2, 5) #rho, tau / no lambda



for (sub in 1:N){
  sub_id <- subject_id[sub]
  
  data_sub <- subset(ra_data, subjID == sub_id)
  
  mle_model_noLA <- optim(param_init, ra_noLA, method="L-BFGS-B", lower=param_noLA_low, upper=param_noLA_up, ra_data = data_sub)
  
  for (i in 1:100) {
    
    param_init <- runif(2)
    
    # Do the MLE again
    temp_noLA <- optim(param_init, ra_noLA, method="L-BFGS-B", lower=param_noLA_low, upper=param_noLA_up, ra_data = data_sub)
    
    # Replace the results if the latest optimization yields better result
    if(temp_noLA$value < mle_model_noLA$value)  mle_model_noLA <- temp_noLA
    
  }
  
  AIC_noLA <- 2 * mle_model_noLA$value + 2 * k_noLA
  BIC_noLA <- 2 * mle_model_noLA$value + k_noLA * log(T) # I wonder we should use T or N / N means a number of subjects, T means a number of trials
  results_noLA$AIC[sub] <- AIC_noLA
  results_noLA$BIC[sub] <- BIC_noLA
}

# results
print(results_noLA)
print(sum(results_noLA$AIC))
print(sum(results_noLA$BIC))

#================================== ra_noRA/// HW 3-2========================
results_noRA <- data.frame(subjID = subject_id, AIC = NA, BIC = NA)

total_noRA <- 0 # for all subjects AIC, BIC
k_noRA <- 2 #every subjects has two parameters


#parameters 
param_init <- runif(2)

#bounds
param_noRA_low <- c(0, 0); param_noRA_up <- c(10, 5) #lambda, tau / no rho



for (sub in 1:N){
  sub_id <- subject_id[sub]
  
  data_sub <- subset(ra_data, subjID == sub_id)
  
  mle_model_noRA <- optim(param_init, ra_noRA, method="L-BFGS-B", lower=param_noRA_low, upper=param_noRA_up, ra_data = data_sub)
  
  for (i in 1:100) {
    
    param_init <- runif(2)
    
    # Do the MLE again
    temp_noRA <- optim(param_init, ra_noRA, method="L-BFGS-B", lower=param_noRA_low, upper=param_noRA_up, ra_data = data_sub)
    
    # Replace the results if the latest optimization yields better result
    if(temp_noRA$value < mle_model_noRA$value)  mle_model_noRA <- temp_noRA
    
  }
  
  AIC_noRA <- 2 * mle_model_noRA$value + 2 * k_noRA
  BIC_noRA <- 2 * mle_model_noRA$value + k_noRA * log(T) # I wonder we should use T or N / N means a number of subjects, T means a number of trials
  results_noRA$AIC[sub] <- AIC_noRA
  results_noRA$BIC[sub] <- BIC_noRA
}

# results
print(results_noRA)
print(sum(results_noRA$AIC))
print(sum(results_noRA$BIC))

#================================== all models /// HW 3-3========================
# we already calculate ra_noLA, ra_noRA, so AIC, BIC of ra_prospect is left
#For consistency, rewrote the prospect model to match the noRA and noLA models

results_prospect <- data.frame(subjID = subject_id, AIC = NA, BIC = NA)

total_prospect <- 0 # for all subjects AIC, BIC
k_prospect <- 3 #every subjects has two parameters


#parameters 
param_init <- runif(3)

#bounds
param_prospect_low <- c(0, 0, 0); param_prospect_up <- c(2, 10, 5) #rho, lambda, tau


for (sub in 1:N){
  sub_id <- subject_id[sub]
  
  data_sub <- subset(ra_data, subjID == sub_id)
  
  mle_model_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub)
  
  for (i in 1:100) {
    
    param_init <- runif(3)
    
    # Do the MLE again
    temp_prospect <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, ra_data = data_sub)
    
    # Replace the results if the latest optimization yields better result
    if(temp_prospect$value < mle_model_prospect$value)  mle_model_prospect <- temp_prospect
    
  }
  
  AIC_prospect <- 2 * mle_model_prospect$value + 2 * k_prospect
  BIC_prospect <- 2 * mle_model_prospect$value + k_prospect * log(T) # I wonder we should use T or N / N means a number of subjects, T means a number of trials
  results_prospect$AIC[sub] <- AIC_prospect
  results_prospect$BIC[sub] <- BIC_prospect
}

# results
print(results_prospect)
print(sum(results_prospect$AIC))
print(sum(results_prospect$BIC))

# comparing three ra_models


all_AIC = round(c(sum(results_prospect$AIC), sum(results_noLA$AIC), sum(results_noRA$AIC)), 3)
all_BIC = round(c(sum(results_prospect$BIC), sum(results_noLA$BIC), sum(results_noRA$BIC)), 3)
names = c("ra_prospect", "ra_noLA", "ra_noRA")

modelcomp_summary = data.frame(Models = names, AIC = all_AIC, BIC = all_BIC)
print(modelcomp_summary)
