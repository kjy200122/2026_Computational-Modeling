#######################################################
## Parm_Estm.R                                     ##
## PSYCH 7695, Spring 2017                           ##
## Maximum Likelihood Estimation (Myung, JMP, 2003)  ##
## By Yun Tang, Psychology, OSU                      ##
##                                                   ##
## Main Program                                      ##
## Code Written on 12/18/2009                        ##
##                                                   ##
## Modified by Joonsuk Park on Jan 28 2015           ##
## Modified by Jay Myung in Feb 2017                 ##
## Modified by Woo-Young Ahn in March 2018           ##
#######################################################

# Loading the (minus) log-likelihood functions
# Please modify the path according to the actual location of the file "MLE_LSE.R"
# e.g., setwd("/Users/youngahn/this-course/")

rm(list=ls())  # clear workspace
graphics.off() # close all figures

set.seed(202604)  # set a seed number for replication
source("MLE_LSE.R")   # source MLE_LSE.R code

##########################
## General Setup        ##
## Data and Parameters  ##
##########################
n_total <- 50 # sample size 
t_int <-  c(0.5, 1,  2,  4,  8, 12, 18, 20) # time interval values
n_corr <- c(44, 35, 27, 24, 17, 15, 13, 10) # number of correct responses
p_corr <- n_corr/n_total # proportion correct

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param1_init <- runif(1)
param2_init <- runif(2)
param3_init <- runif(3)

param_pow1_low <- 0;         param_pow1_up <- 5;
param_pow2_low <- c(0, 0);   param_pow2_up <- c(1, 5);
param_exp1_low <- 0;         param_exp1_up <- 5;
param_exp2_low <- c(0, 0);   param_exp2_up <- c(1, 5);
param_expow_low <- c(0,0,0); param_expow_up <- c(1,5,5);
param_hyp1_low <- 0;         param_hyp1_up <- 1;
param_hyp2_low <- c(0,0);    param_hyp2_up <- c(1,1);

##########################
## MLE                  ##
##########################

# Custom function for optimization with default settings
optimize <- function(init, func, lower, upper, ...) {
  optim(init, func, lower = lower, upper = upper, method = 'L-BFGS-B',
        int = t_int, n = n_total, x = n_corr, ...)
}
# Call general purpose optimization rountine
mle_model_pow1 <- optimize(param1_init, mle_pow1, param_pow1_low, param_pow1_up)
mle_model_pow2 <- optimize(param2_init, mle_pow2, param_pow2_low, param_pow2_up)
mle_model_exp1 <- optimize(param1_init, mle_exp1, param_exp1_low, param_exp1_up)
mle_model_exp2 <- optimize(param2_init, mle_exp2, param_exp2_low, param_exp2_up)
mle_model_expow <- optimize(param3_init, mle_expow, param_expow_low, param_expow_up)
mle_model_hyp1 <- optimize(param1_init, mle_hyp1, param_hyp1_low, param_hyp1_up)
mle_model_hyp2 <- optimize(param2_init, mle_hyp2, param_hyp2_low, param_hyp2_up)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
  temp_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
  temp_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
  temp_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
  temp_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp2 <- optim(param2_init, mle_hyp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)
  
  # Replace the results if the latest optimization yields better result
  if(temp_pow1$value < mle_model_pow1$value) mle_model_pow1 <- temp_pow1 
  if(temp_pow2$value < mle_model_pow2$value) mle_model_pow2 <- temp_pow2
  if(temp_exp1$value < mle_model_exp1$value) mle_model_exp1 <- temp_exp1
  if(temp_exp2$value < mle_model_exp2$value) mle_model_exp2 <- temp_exp2
  if(temp_expow$value < mle_model_expow$value) mle_model_expow <- temp_expow
  if(temp_hyp1$value < mle_model_hyp1$value) mle_model_hyp1 <- temp_hyp1
  if(temp_hyp2$value < mle_model_hyp2$value) mle_model_hyp2 <- temp_hyp2
}

# Save the MLE parameter estimates
parm_pow1 <- mle_model_pow1$par
parm_pow2 <- mle_model_pow2$par
parm_exp1 <- mle_model_exp1$par
parm_exp2 <- mle_model_exp2$par
parm_expow <- mle_model_expow$par
parm_hyp1 <- mle_model_hyp1$par
parm_hyp2 <- mle_model_hyp2$par

# MLE predictions
p_prd_pow1  <- compute_p_pow1(t_int, parm_pow1)
p_prd_pow2  <- compute_p_pow2(t_int, parm_pow2)
p_prd_exp1  <- compute_p_exp1(t_int, parm_exp1)
p_prd_exp2  <- compute_p_exp2(t_int, parm_exp2)
p_prd_expow <- compute_p_expow(t_int, parm_expow)
p_prd_hyp1  <- compute_p_hyp1(t_int, parm_hyp1)
p_prd_hyp2  <- compute_p_hyp2(t_int, parm_hyp2)

# Proportion of the explained variances for each model
calculte_r2 <- function(p_corr, p_prd){
  r2 = 1-sum((p_corr-p_prd)^2)/sum((p_corr-mean(p_corr))^2)
  return(r2)
}
r2_pow1 = calculte_r2(p_corr,p_prd_pow1)
r2_pow2 = calculte_r2(p_corr,p_prd_pow2)
r2_exp1 = calculte_r2(p_corr,p_prd_exp1)
r2_exp2 = calculte_r2(p_corr,p_prd_exp2)
r2_expow = calculte_r2(p_corr,p_prd_expow)
r2_hyp1 = calculte_r2(p_corr,p_prd_hyp1)
r2_hyp2 = calculte_r2(p_corr,p_prd_hyp2)

# Generate MLE summary
model_names = c("POW1","POW2","EXP1", "EXP2","EXPOW","HYP1","HYP2")

minus_loglik_mle = round(c(mle_model_pow1$value,
                           mle_model_pow2$value,
                           mle_model_exp1$value,
                           mle_model_exp2$value,
                           mle_model_expow$value,
                           mle_model_hyp1$value,
                           mle_model_hyp2$value), 3)

r2_mle <- round(c(r2_pow1,
                  r2_pow2,
                  r2_exp1,
                  r2_exp2,
                  r2_expow,
                  r2_hyp1,
                  r2_hyp2), 3)

mle_summary = data.frame(Models = model_names,
                         loglik = -minus_loglik_mle,
                         r2 = r2_mle)

pars_mle <- round(cbind(mle_model_pow1$par,
                        mle_model_pow2$par,
                        mle_model_exp1$par,
                        mle_model_exp2$par,
                        mle_model_expow$par,
                        mle_model_hyp1$par,
                        mle_model_hyp2$par), 3)

dimnames(pars_mle) = list(c('par1', 'par2', 'par3'), model_names)


##########################
## LSE                  ##
##########################

# Call general purpose optimization rountine
lse_model_pow1 <- optimize(param1_init, lse_pow1, param_pow1_low, param_pow1_up)
lse_model_pow2 <- optimize(param2_init, lse_pow2, param_pow2_low, param_pow2_up)
lse_model_exp1 <- optimize(param1_init, lse_exp1, param_exp1_low, param_exp1_up)
lse_model_exp2 <- optimize(param2_init, lse_exp2, param_exp2_low, param_exp2_up)
lse_model_expow <- optimize(param3_init, lse_expow, param_expow_low, param_expow_up)
lse_model_hyp1 <- optimize(param1_init, lse_hyp1, param_hyp1_low, param_hyp1_up)
lse_model_hyp2 <- optimize(param2_init, lse_hyp2, param_hyp2_low, param_hyp2_up)


# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the LSE again
  temp_pow1  <- optimize(param1_init, lse_pow1,  param_pow1_low,  param_pow1_up)
  temp_pow2  <- optimize(param2_init, lse_pow2,  param_pow2_low,  param_pow2_up)
  temp_exp1  <- optimize(param1_init, lse_exp1,  param_exp1_low,  param_exp1_up)
  temp_exp2  <- optimize(param2_init, lse_exp2,  param_exp2_low,  param_exp2_up)
  temp_expow <- optimize(param3_init, lse_expow, param_expow_low, param_expow_up)
  temp_hyp1  <- optimize(param1_init, lse_hyp1,  param_hyp1_low,  param_hyp1_up)
  temp_hyp2  <- optimize(param2_init, lse_hyp2,  param_hyp2_low,  param_hyp2_up)
  
  # Replace the results if the latest optimization yields better result
  if (temp_pow1$value  < lse_model_pow1$value)  lse_model_pow1  <- temp_pow1
  if (temp_pow2$value  < lse_model_pow2$value)  lse_model_pow2  <- temp_pow2
  if (temp_exp1$value  < lse_model_exp1$value)  lse_model_exp1  <- temp_exp1
  if (temp_exp2$value  < lse_model_exp2$value)  lse_model_exp2  <- temp_exp2
  if (temp_expow$value < lse_model_expow$value) lse_model_expow <- temp_expow
  if (temp_hyp1$value  < lse_model_hyp1$value)  lse_model_hyp1  <- temp_hyp1
  if (temp_hyp2$value  < lse_model_hyp2$value)  lse_model_hyp2  <- temp_hyp2
}

# Save the LSE parameter estimates
parm_lse_pow1  <- lse_model_pow1$par
parm_lse_pow2  <- lse_model_pow2$par
parm_lse_exp1  <- lse_model_exp1$par
parm_lse_exp2  <- lse_model_exp2$par
parm_lse_expow <- lse_model_expow$par
parm_lse_hyp1  <- lse_model_hyp1$par
parm_lse_hyp2  <- lse_model_hyp2$par

# lse predictions
p_prd_pow1  <- compute_p_pow1(t_int, parm_lse_pow1)
p_prd_pow2  <- compute_p_pow2(t_int, parm_lse_pow2)
p_prd_exp1  <- compute_p_exp1(t_int, parm_lse_exp1)
p_prd_exp2  <- compute_p_exp2(t_int, parm_lse_exp2)
p_prd_expow <- compute_p_expow(t_int, parm_lse_expow)
p_prd_hyp1  <- compute_p_hyp1(t_int, parm_lse_hyp1)
p_prd_hyp2  <- compute_p_hyp2(t_int, parm_lse_hyp2)

# Proportion of the explained variances for each model
calculte_r2 <- function(p_corr, p_prd){
  r2 = 1-sum((p_corr-p_prd)^2)/sum((p_corr-mean(p_corr))^2)
  return(r2)
}
r2_pow1 = calculte_r2(p_corr,p_prd_pow1)
r2_pow2 = calculte_r2(p_corr,p_prd_pow2)
r2_exp1 = calculte_r2(p_corr,p_prd_exp1)
r2_exp2 = calculte_r2(p_corr,p_prd_exp2)
r2_expow = calculte_r2(p_corr,p_prd_expow)
r2_hyp1 = calculte_r2(p_corr,p_prd_hyp1)
r2_hyp2 = calculte_r2(p_corr,p_prd_hyp2)

# Generate summary
model_names <- c('POW1', 'POW2', 'EXP1', 'EXP2', 'EXPOW', 'HYP1', 'HYP2')

sse = round(c(lse_model_pow1$value,
                           lse_model_pow2$value,
                           lse_model_exp1$value,
                           lse_model_exp2$value,
                           lse_model_expow$value,
                           lse_model_hyp1$value,
                           lse_model_hyp2$value), 3)

r2_lse <- round(c(r2_pow1,
                  r2_pow2,
                  r2_exp1,
                  r2_exp2,
                  r2_expow,
                  r2_hyp1,
                  r2_hyp2), 3)

lse_summary = data.frame(Models = model_names,
                         SSE = sse,
                         r2 = r2_lse)

pars_lse <- round(cbind(lse_model_pow1$par,
                        lse_model_pow2$par,
                        lse_model_exp1$par,
                        lse_model_exp2$par,
                        lse_model_expow$par,
                        lse_model_hyp1$par,
                        lse_model_hyp2$par), 3)

dimnames(pars_lse) = list(c('par1', 'par2', 'par3'), model_names)

# Load libraries for visualization
library(ggplot2)
library(tidyverse)

# Visualization of MLE results
x <- seq(0, 20, 0.05)

# Compute predictions for each model
p_pow1 <- (x + 1) ^ (-parm_pow1[1])
p_pow2 <- parm_pow2[1] * (x + 1) ^ (-parm_pow2[2])
p_exp1 <- exp(-parm_exp1[1] * x)
p_exp2 <- parm_exp2[1] * exp(-parm_exp2[2] * x)
p_expow <- parm_expow[1] * exp(-parm_expow[2] * x) * (x + 1) ^ (-parm_expow[3])
p_hyp1 <- 1 / (1 + parm_hyp1[1] * x)
p_hyp2 <- parm_hyp2[1] / (1 + parm_hyp2[2] * x)

# Combine into a long-format data frame
simulation_data <- data.frame(x, p_pow1, p_pow2, p_exp1, p_exp2, p_expow, p_hyp1, p_hyp2) %>%
  pivot_longer(cols = starts_with("p"), names_to = "Model", values_to = "Pred_Pr")

# Add linetypes and model labels
simulation_data$linetype <- as.factor(rep(c(1, 2, 1, 2, 3, 1, 2), each = length(x)))
simulation_data$Model <- factor(simulation_data$Model, 
                                levels = c("p_pow1", "p_pow2", "p_exp1", "p_exp2", "p_expow", "p_hyp1", "p_hyp2"))
cols <- c("red", "orange", "yellow", "green", "blue", "navy", "purple")
text <- c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")

# Plot MLE results
plot_mle <- ggplot(simulation_data, aes(x = x)) +
  geom_line(aes(y = Pred_Pr, colour = Model, linetype = linetype)) +
  geom_point(data = data.frame(t_int, p_corr), aes(x = t_int, y = p_corr)) +
  labs(title = "Probability of recall predicted by each model (MLE)", 
       x = "Time (t)", 
       y = "Proportion Correct") +
  scale_colour_manual(name = "Models", values = cols, labels = text) +
  scale_linetype(guide = "none")

print(plot_mle)

# Visualization of LSE results

# Compute predictions for each model with LSE parameters
p_lse_pow1 <- (x + 1) ^ (-parm_lse_pow1[1])
p_lse_pow2 <- parm_lse_pow2[1] * (x + 1) ^ (-parm_lse_pow2[2])
p_lse_exp1 <- exp(-parm_lse_exp1[1] * x)
p_lse_exp2 <- parm_lse_exp2[1] * exp(-parm_lse_exp2[2] * x)
p_lse_expow <- parm_lse_expow[1] * exp(-parm_lse_expow[2] * x) * (x + 1) ^ (-parm_lse_expow[3])
p_lse_hyp1 <- 1 / (1 + parm_lse_hyp1[1] * x)
p_lse_hyp2 <- parm_lse_hyp2[1] / (1 + parm_lse_hyp2[2] * x)

# Combine into a long-format data frame
lse_simulation_data <- data.frame(x, p_lse_pow1, p_lse_pow2, p_lse_exp1, p_lse_exp2, p_lse_expow, p_lse_hyp1, p_lse_hyp2) %>%
  pivot_longer(cols = starts_with("p"), names_to = "Model", values_to = "Pred_Pr")

# Add linetypes and model labels
lse_simulation_data$linetype <- as.factor(rep(c(1, 2, 1, 2, 3, 1, 2), each = length(x)))
lse_simulation_data$Model <- factor(lse_simulation_data$Model, 
                                    levels = c("p_lse_pow1", "p_lse_pow2", "p_lse_exp1", "p_lse_exp2", "p_lse_expow", "p_lse_hyp1", "p_lse_hyp2"))

# Plot LSE results
plot_lse <- ggplot(lse_simulation_data, aes(x = x)) +
  geom_line(aes(y = Pred_Pr, colour = Model, linetype = linetype)) +
  geom_point(data = data.frame(t_int, p_corr), aes(x = t_int, y = p_corr)) +
  labs(title = "Probability of recall predicted by each model (LSE)", 
       x = "Time (t)", 
       y = "Proportion Correct") +
  scale_colour_manual(name = "Models", values = cols, labels = text) +
  scale_linetype(guide = "none")

print(plot_lse)
