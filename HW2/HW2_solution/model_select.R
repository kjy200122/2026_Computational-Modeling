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

set.seed(202604)  # set a seed number for replication
source("MLE.R")   # source MLE.R code

##########################
## General Setup        ##
## Data and Parameters  ##
##########################
n_total <- 50 # sample size
t_int <-  c(0.5, 1,  2,  4,  8, 12, 18, 20) # time interval values
n_corr <- c(44, 35, 27, 24, 17, 15, 13, 11) # number of correct responses
p_corr <- n_corr / n_total # proportion correct

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param1_init <- runif(1)
param2_init <- runif(2)
param3_init <- runif(3)

param_pow1_low  <- c(0);       param_pow1_up  <- c(5);
param_pow2_low  <- c(0, 0);    param_pow2_up  <- c(1, 5);
param_exp1_low  <- c(0);       param_exp1_up  <- c(5);
param_exp2_low  <- c(0, 0);    param_exp2_up  <- c(1, 5);
param_expow_low <- c(0, 0, 0); param_expow_up <- c(1, 5, 5);
param_hyp1_low  <- c(0);       param_hyp1_up  <- c(1);
param_hyp2_low  <- c(0, 0);    param_hyp2_up  <- c(1, 1);

# Custom function for optimization with default settings
optimize_hw2 <- function(init, func, lower, upper, ...) {
  optim(init, func, lower = lower, upper = upper, method = 'L-BFGS-B',
        int = t_int, n = n_total, x = n_corr, ...)
}

##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_pow1 <- optimize_hw2(param1_init, mle_pow1, param_pow1_low, param_pow1_up)
mle_model_pow2 <- optimize_hw2(param2_init, mle_pow2, param_pow2_low, param_pow2_up)
mle_model_exp1 <- optimize_hw2(param1_init, mle_exp1, param_exp1_low, param_exp1_up)
mle_model_exp2 <- optimize_hw2(param2_init, mle_exp2, param_exp2_low, param_exp2_up)
mle_model_expow <- optimize_hw2(param3_init, mle_expow, param_expow_low, param_expow_up)
mle_model_hyp1 <- optimize_hw2(param1_init, mle_hyp1, param_hyp1_low, param_hyp1_up)
mle_model_hyp2 <- optimize_hw2(param2_init, mle_hyp2, param_hyp2_low, param_hyp2_up)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_pow1  <- optimize_hw2(param1_init, mle_pow1,  param_pow1_low,  param_pow1_up)
  temp_pow2  <- optimize_hw2(param2_init, mle_pow2,  param_pow2_low,  param_pow2_up)
  temp_exp1  <- optimize_hw2(param1_init, mle_exp1,  param_exp1_low,  param_exp1_up)
  temp_exp2  <- optimize_hw2(param2_init, mle_exp2,  param_exp2_low,  param_exp2_up)
  temp_expow <- optimize_hw2(param3_init, mle_expow, param_expow_low, param_expow_up)
  temp_hyp1  <- optimize_hw2(param1_init, mle_hyp1,  param_hyp1_low,  param_hyp1_up)
  temp_hyp2  <- optimize_hw2(param2_init, mle_hyp2,  param_hyp2_low,  param_hyp2_up)  
  
  # Replace the results if the latest optimization yields better result
  if (temp_pow1$value  < mle_model_pow1$value)  mle_model_pow1  <- temp_pow1
  if (temp_pow2$value  < mle_model_pow2$value)  mle_model_pow2  <- temp_pow2
  if (temp_exp1$value  < mle_model_exp1$value)  mle_model_exp1  <- temp_exp1
  if (temp_exp2$value  < mle_model_exp2$value)  mle_model_exp2  <- temp_exp2
  if (temp_expow$value < mle_model_expow$value) mle_model_expow <- temp_expow
  if (temp_hyp1$value  < mle_model_hyp1$value)  mle_model_hyp1  <- temp_hyp1
  if (temp_hyp2$value  < mle_model_hyp2$value)  mle_model_hyp2  <- temp_hyp2
}

# Save the MLE parameter estimates
parm_pow1  <- mle_model_pow1$par
parm_pow2  <- mle_model_pow2$par
parm_exp1  <- mle_model_exp1$par
parm_exp2  <- mle_model_exp2$par
parm_expow <- mle_model_expow$par
parm_hyp1  <- mle_model_hyp1$par
parm_hyp2  <- mle_model_hyp2$par

# number of parameters (k_modelName)
k_pow1 <- 1
k_pow2 <- 2
k_exp1 <- 1
k_exp2 <- 2
k_expow <- 3
k_hyp1 <- 1
k_hyp2 <- 2

# number of data points (N)
N = length(p_corr) 

# Compute AIC = -2*log(lik) + 2*K
compute_AIC <- function(log_lik, k){
  AIC <- 2*log_lik + 2*k
  return(AIC)
}

AIC_pow1 <- compute_AIC(mle_model_pow1$value, k_pow1) # mle_model_pow2$value --> -log(lik)
AIC_pow2 <- compute_AIC(mle_model_pow2$value, k_pow2)
AIC_exp1 <- compute_AIC(mle_model_exp1$value, k_exp1)
AIC_exp2 <- compute_AIC(mle_model_exp2$value, k_exp2)
AIC_expow <- compute_AIC(mle_model_expow$value, k_expow)
AIC_hyp1 <- compute_AIC(mle_model_hyp1$value, k_hyp1)
AIC_hyp2 <- compute_AIC(mle_model_hyp2$value, k_hyp2)

# Compute BIC = -2*log(lik) + K*log(N)
compute_BIC <- function(log_lik, k, N){
  BIC <- 2*log_lik + k*log(N)
  return(BIC)
}

BIC_pow1 <- compute_BIC(mle_model_pow1$value, k_pow1, N) # mle_model_pow2$value --> -log(lik)
BIC_pow2 <- compute_BIC(mle_model_pow2$value, k_pow2, N)
BIC_exp1 <- compute_BIC(mle_model_exp1$value, k_exp1, N)
BIC_exp2 <- compute_BIC(mle_model_exp2$value, k_exp2, N)
BIC_expow <- compute_BIC(mle_model_expow$value, k_expow, N)
BIC_hyp1 <- compute_BIC(mle_model_hyp1$value, k_hyp1, N)
BIC_hyp2 <- compute_BIC(mle_model_hyp2$value, k_hyp2, N)

# Generate summary
model_names <- c('POW1', 'POW2', 'EXP1', 'EXP2', 'EXPOW', 'HYP1', 'HYP2')

all_AIC = round(c(AIC_pow1, 
                  AIC_pow2,
                  AIC_exp1,
                  AIC_exp2,
                  AIC_expow,
                  AIC_hyp1,
                  AIC_hyp2), 3)

all_BIC = round(c(BIC_pow1, 
                  BIC_pow2,
                  BIC_exp1,
                  BIC_exp2,
                  BIC_expow,
                  BIC_hyp1,
                  BIC_hyp2), 3)

modelcomp_summary = data.frame(Models = model_names,
                               AIC = all_AIC,
                               BIC = all_BIC
                               )
print(modelcomp_summary)
