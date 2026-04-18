################################################################################
## model_select.R                                                             ##
## PSYCH 7695, Spring 2017                                                    ##
## Maximum Likelihood Estimation (Myung, JMP, 2003)                           ##
## By Yun Tang, Psychology, OSU                                               ##
##                                                                            ##
## Main Program                                                               ##
## Code Written on 12/18/2009                                                 ##
##                                                                            ##
## Modified by Joonsuk Park on Jan 28 2015                                    ##
## Modified by Jay Myung in Feb 2017                                          ##
## Modified by Woo-Young Ahn in March 2018                                    ##
## Modified by Jaeyeong Yang in March 2019                                    ##
## Modified by Woo-Young Ahn in May 2021                                      ##
################################################################################

# Loading the (minus) log-likelihood functions
rm(list=ls())  # clear workspace
graphics.off() # close all figures

setwd('~/Downloads/Homework1_소은솔')

set.seed(202604)  # set a seed number for replication
source("MLE_LSE.R")   # source MLE._LSE.R code you programmed for HW1

##########################
## General Setup        ##
## Data and Parameters  ##
##########################
n_total <- 50 # sample size
t_int <-  c(0.5, 1,  2,  4,  8, 12, 18, 20) # time interval values
n_corr <- c(44, 35, 27, 24, 17, 15, 13, 11) # number of correct responses
p_corr <- n_corr/n_total # proportion correct

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param1_init <- runif(1)
param2_init <- runif(2)
param3_init <- runif(3)

param_pow2_low <- c(0, 0); param_pow2_up <- c(1, 5);  # lower and upper bounds of POW2 model (0<a<1, 0<b<5)
param_exp2_low <- c(0, 0); param_exp2_up <- c(1, 5);  # lower and upper bounds of EXP2 model (0<a<1, 0<b<5)
param_pow1_low <- c(0); param_pow1_up <- c(5);
param_exp1_low <- c(0); param_exp1_up <- c(5);
param_expow_low <- c(0, 0, 0); param_expow_up <- c(1, 5, 5);
param_hyp1_low <- c(0); param_hyp1_up <- c(1);
param_hyp2_low <- c(0,0); param_hyp2_up <- c(1,1);

##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
mle_model_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp2 <- optim(param2_init, mle_hyp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
  temp_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
  temp_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
  temp_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
  temp_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp2 <- optim(param2_init, mle_hyp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)
  
  # Replace the results if the latest optimization yields better result
  if(temp_pow2$value < mle_model_pow2$value) mle_model_pow2 <- temp_pow2  
  if(temp_exp2$value < mle_model_exp2$value) mle_model_exp2 <- temp_exp2
  if(temp_pow1$value < mle_model_pow1$value) mle_model_pow1 <- temp_pow1
  if(temp_exp1$value < mle_model_exp1$value) mle_model_exp1 <- temp_exp1
  if(temp_expow$value < mle_model_expow$value) mle_model_expow <- temp_expow
  if(temp_hyp1$value < mle_model_hyp1$value) mle_model_hyp1 <- temp_hyp1
  if(temp_hyp2$value < mle_model_hyp2$value) mle_model_hyp2 <- temp_hyp2
}

# Save the MLE parameter estimates
parm_pow2 <- mle_model_pow2$par
parm_exp2 <- mle_model_exp2$par
parm_pow1 <- mle_model_pow1$par
parm_exp1 <- mle_model_exp1$par
parm_expow <- mle_model_expow$par
parm_hyp1 <- mle_model_hyp1$par
parm_hyp2 <- mle_model_hyp2$par

# number of parameters (k_modelName)
k_pow2 <- 2
k_exp2 <- 2
k_pow1 <- 1
k_exp1 <- 1
k_expow <- 3
k_hyp1 <- 1
k_hyp2 <- 2

# number of data points (N)
N = length(p_corr) 

# Compute AIC = -2*log(lik) + 2*K
AIC_pow2 <- 2*mle_model_pow2$value + 2*k_pow2   # mle_model_pow2$value --> -log(lik)
AIC_exp2 <- 2*mle_model_exp2$value + 2*k_exp2
AIC_pow1 <- 2*mle_model_pow1$value + 2*k_pow1
AIC_exp1 <- 2*mle_model_exp1$value + 2*k_exp1
AIC_expow <- 2*mle_model_expow$value + 2*k_expow
AIC_hyp1 <- 2*mle_model_hyp1$value + 2*k_hyp1
AIC_hyp2 <- 2*mle_model_hyp2$value + 2*k_hyp2

# Compute BIC = -2*log(lik) + K*log(N)
BIC_pow2 <- 2*mle_model_pow2$value + k_pow2*log(N)   # mle_model_pow2$value --> -log(lik)
BIC_exp2 <- 2*mle_model_exp2$value + k_exp2*log(N)
BIC_pow1 <- 2*mle_model_pow1$value + k_pow1*log(N)
BIC_exp1 <- 2*mle_model_exp1$value + k_exp1*log(N)
BIC_expow <- 2*mle_model_expow$value + k_expow*log(N)
BIC_hyp1 <- 2*mle_model_hyp1$value + k_hyp1*log(N)
BIC_hyp2 <- 2*mle_model_hyp2$value + k_hyp2*log(N)

# Generate summary
all_AIC = round(c(AIC_pow2, AIC_exp2, AIC_pow1, AIC_exp1, AIC_expow, AIC_hyp1, AIC_hyp2), 3)
all_BIC = round(c(BIC_pow2, BIC_exp2, BIC_pow1, BIC_exp1, BIC_expow, BIC_hyp1, BIC_hyp2), 3)
names = c("POW2", "EXP2", "POW1", "EXP1","EXPOW", "HYP1", "HYP2")

modelcomp_summary = data.frame(Models = names, AIC = all_AIC, BIC = all_BIC)

print(modelcomp_summary)

########################################################
#Question 2
#for subject 2
########################################################

rm(list=ls())  # clear workspace
graphics.off() # close all figures

setwd('~/Downloads/HW2')
data <- read.table("ra_exampleData.txt", header = TRUE)

set.seed(202604)

#filter data for subj 2
subj2_data <- subset(data, subjID == 2)

#MLE function 
ra_prospect <- function(params, cert, gain, loss, gamble){
  rho<- params[1]
  lambda <- params[2]
  tau <- params[3]
  
  N = 5  # number of subjects
  T = 140 # number of trials per subject
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize. 
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)}

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param3_init <- runif(3)
#lower and upper bounds 
param_prospect_low <- c(0, 0, 0); param_prospect_up <- c(2, 10, 5)

#general purpose optimisation routine 
mle_model_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                            cert=subj2_data$cert, gain=subj2_data$gain, loss=subj2_data$loss, gamble=subj2_data$gamble)

for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param3_init <- runif(3); 
  
  # Do the MLE again
  temp_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                         cert=subj2_data$cert, gain=subj2_data$gain, loss=subj2_data$loss, gamble=subj2_data$gamble)
  
  # Replace the results if the latest optimization yields better result
  if(temp_prospect$value < mle_model_prospect$value) mle_model_prospect <- temp_prospect  
}

# Save the MLE parameter estimates
parm_final <- round(mle_model_prospect$par, 3)
mle_summary = data.frame(
  Subject = 2,
  rho     = parm_final[1],
  lambda  = parm_final[2],
  tau     = parm_final[3]
)

print('- MLE results for Prospect Model (Subj 2) ------------')
print(mle_summary)

###################################################
#For all subjects
##################################################

all_subjects <- c(2,3,4,6,7)
all_results <- data.frame()

for (s in all_subjects) {
  
  subj_data <- subset(data, subjID == s)
  param3_init<-runif(3)
  mle_model_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                              cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
  for (i in 1:100) {
    # Re-generate random inits. Is it the best way to do this?
    param3_init <- runif(3); 
    
    # Do the MLE again
    temp_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                           cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_prospect$value < mle_model_prospect$value) mle_model_prospect <- temp_prospect  
  }
  final <- round(mle_model_prospect$par, 3)
  mle_summary = data.frame(
    Subject = s,
    rho     = final[1],
    lambda  = final[2],
    tau     = final[3]
  )
  all_results <-rbind(all_results, mle_summary)
}

print('-Results for all subjects-----------------------')
print(all_results)

##############################################################
#Question 3 (optional)
##############################################################

set.seed(202604)

#function without loss aversion parameter 
ra_noLA <- function(params, cert, gain, loss, gamble){
  rho<- params[1]
  tau <- params[2]
  
  N = 5  # number of subjects
  T = 140 # number of trials per subject
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize. 
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)}

#############################################################
#AIC and BIC for all participants 
#############################################################
model_comparison_results <- data.frame()

for (s in all_subjects) {
  
  subj_data <- subset(data, subjID == s)
  param2_init<-runif(2)
  
  param_noLA_low <- c(0, 0); param_noLA_up <- c(2, 5)
  mle_model_noLA <- optim(param2_init, ra_noLA, method="L-BFGS-B", lower=param_noLA_low, upper=param_noLA_up, 
                              cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
  for (i in 1:100) {
    # Re-generate random inits. Is it the best way to do this?
    param2_init <- runif(2); 
    
    # Do the MLE again
    temp_noLA <- optim(param2_init, ra_noLA, method="L-BFGS-B", lower=param_noLA_low, upper=param_noLA_up, 
                           cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_noLA$value < mle_model_noLA$value) mle_model_noLA <- temp_noLA  
  }
  
  k_noLA <- 2 
  datpoints = 140
  AIC_noLA <- 2*mle_model_noLA$value + 2*k_noLA
  BIC_noLA <- 2*mle_model_noLA$value + k_noLA*log(datpoints)
  
  subj_metrics <- data.frame(
    SubjectID = s,
    Model     = "noLA",
    AIC       = round(AIC_noLA, 3),
    BIC       = round(BIC_noLA, 3)
  )
  model_comparison_results <- rbind(model_comparison_results, subj_metrics)
}

total_row <- data.frame(
  SubjectID = "Total",
  Model = "noLA",
  AIC = round(sum(model_comparison_results$AIC), 3),
  BIC = round(sum(model_comparison_results$BIC), 3)
)

final_summary_table <- rbind(model_comparison_results, total_row)

print("--- Model Comparison Summary (noLA Model) ---")
print(final_summary_table, row.names = FALSE)

########################################################
#function without risk preference parameter 
########################################################

ra_noRA <- function(params, cert, gain, loss, gamble){
  lambda<- params[1]
  tau <- params[2]
  
  N = 5  # number of subjects
  T = 140 # number of trials per subject
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize. 
  for (t in 1:T) {
    evSafe   = cert[t]
    evGamble = 0.5*(gain[t] - lambda*abs(loss[t])) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)}

####################AIC and BIC results ######################
model_comparison_results <- data.frame()

for (s in all_subjects) {
  
  subj_data <- subset(data, subjID == s)
  param2_init<-runif(2)
  
  param_noRA_low <- c(0, 0); param_noRA_up <- c(10, 5)
  mle_model_noRA <- optim(param2_init, ra_noRA, method="L-BFGS-B", lower=param_noRA_low, upper=param_noRA_up, 
                          cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
  for (i in 1:100) {
    # Re-generate random inits. Is it the best way to do this?
    param2_init <- runif(2); 
    
    # Do the MLE again
    temp_noRA <- optim(param2_init, ra_noRA, method="L-BFGS-B", lower=param_noRA_low, upper=param_noRA_up, 
                       cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_noRA$value < mle_model_noRA$value) mle_model_noRA <- temp_noRA  
  }
  
  k_noRA <- 2 
  datpoints = 140
  AIC_noRA <- 2*mle_model_noRA$value + 2*k_noRA
  BIC_noRA <- 2*mle_model_noRA$value + k_noRA*log(datpoints)
  
  subj_metrics <- data.frame(
    SubjectID = s,
    Model     = "noRA",
    AIC       = round(AIC_noRA, 3),
    BIC       = round(BIC_noRA, 3)
  )
  model_comparison_results <- rbind(model_comparison_results, subj_metrics)
}

total_row <- data.frame(
  SubjectID = "Total",
  Model = "noRA",
  AIC = round(sum(model_comparison_results$AIC), 3),
  BIC = round(sum(model_comparison_results$BIC), 3)
)

final_summary_table <- rbind(model_comparison_results, total_row)

print("--- Model Comparison Summary (noRA Model) ---")
print(final_summary_table, row.names = FALSE)

###############################AIC and BIC for all three parameters#################

model_comparison_results <- data.frame()

for (s in all_subjects) {
  
  subj_data <- subset(data, subjID == s)
  param3_init<-runif(3)
  
  param_prospect_low <- c(0, 0, 0); param_prospect_up <- c(2, 10, 5)
  mle_model_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                          cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
  for (i in 1:100) {
    # Re-generate random inits. Is it the best way to do this?
    param3_init <- runif(3); 
    
    # Do the MLE again
    temp_prospect <- optim(param3_init, ra_prospect, method="L-BFGS-B", lower=param_prospect_low, upper=param_prospect_up, 
                       cert=subj_data$cert, gain=subj_data$gain, loss=subj_data$loss, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result
    if(temp_prospect$value < mle_model_prospect$value) mle_model_prospect <- temp_prospect  
  }
  
  k_prospect <- 3 
  datpoints = 140
  AIC_prospect <- 2*mle_model_prospect$value + 2*k_prospect
  BIC_prospect <- 2*mle_model_prospect$value + k_prospect*log(datpoints)
  
  subj_metrics <- data.frame(
    SubjectID = s,
    Model     = "ra_prospect",
    AIC       = round(AIC_prospect, 3),
    BIC       = round(BIC_prospect, 3)
  )
  model_comparison_results <- rbind(model_comparison_results, subj_metrics)
}

total_row <- data.frame(
  SubjectID = "Total",
  Model = "ra_prospect",
  AIC = round(sum(model_comparison_results$AIC), 3),
  BIC = round(sum(model_comparison_results$BIC), 3)
)

final_summary_table <- rbind(model_comparison_results, total_row)

print("--- Model Comparison Summary (Prospect Model) ---")
print(final_summary_table, row.names = FALSE)




