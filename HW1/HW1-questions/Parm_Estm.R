#######################################################
## Parm_Estm.R                                       ##
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

param_pow2_low <- c(0, 0); param_pow2_up <- c(1, 5);  # lower and upper bounds of POW2 model (0<a<1, 0<b<3)
param_exp2_low <- c(0, 0); param_exp2_up <- c(1, 5);  # lower and upper bounds of EXP2 model (0<a<1, 0<b<3)

##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
  temp_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
  
  # Replace the results if the latest optimization yields better result
  if(temp_pow2$value < mle_model_pow2$value) mle_model_pow2 <- temp_pow2  
  if(temp_exp2$value < mle_model_exp2$value) mle_model_exp2 <- temp_exp2
}

# Save the MLE parameter estimates
parm_pow2 <- mle_model_pow2$par
parm_exp2 <- mle_model_exp2$par

# MLE predictions
p_prd_pow2 <- parm_pow2[1]*(t_int+1)^(-parm_pow2[2])
p_prd_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*t_int)

# Proportion of the explained variances for each model
r2_pow2 = 1-sum((p_corr-p_prd_pow2)^2)/sum((p_corr-mean(p_corr))^2)
r2_exp2 = 1-sum((p_corr-p_prd_exp2)^2)/sum((p_corr-mean(p_corr))^2)

# Generate summary
minus_loglik_MLE = round(c(mle_model_pow2$value, mle_model_exp2$value), 3)
r2_mle <- round(c(r2_pow2, r2_exp2), 3)  # round off the value, up to 3 digits
names = c("POW2", "EXP2")
pars_mle <- round(cbind(mle_model_pow2$par, mle_model_exp2$par),3)
dimnames(pars_mle) = list(c('par1', 'par2'),c('POW2', 'EXP2'))

mle_summary = data.frame(Models = names, loglik = - minus_loglik_MLE, r2 = r2_mle)

# Plot the MLE results using simple R graphics. Use ggplots for fancy graphs
x <- seq(0,20, 0.05)
p_pow2 <- parm_pow2[1]*(x+1)^(-parm_pow2[2])
p_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*x)
plot(x, p_pow2, ylim=c(0,1), xlab='Time t', ylab='Proportion Correct', 
     main="MLE results", type='l', lwd=4)
lines(x, p_exp2, lwd=2, lty='dashed', col='red')
points(t_int, p_corr, pch=19, cex=1.5)
text = c("POW2", "EXP2")
ltys = c('solid','dashed')
cols = c('black','red')
legend(14, 1, text, lwd=2, cex=1, lty=ltys, col=cols)

# print maximized likehood values
print('- MLE results ------------')
print(mle_summary,4)

# print bet-fit parameter values
print('- Best-fit parameters --------')
print(pars_mle,4)
