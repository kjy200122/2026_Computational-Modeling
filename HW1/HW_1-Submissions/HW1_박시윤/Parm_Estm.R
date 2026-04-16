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

param_pow1_low <- c(0); param_pow1_up <- c(5);
param_pow2_low <- c(0, 0); param_pow2_up <- c(1, 5);  # lower and upper bounds of POW2 model (0<a<1, 0<b<3)
param_exp1_low <- c(0); param_exp1_up <- c(5);
param_exp2_low <- c(0, 0); param_exp2_up <- c(1, 5);  # lower and upper bounds of EXP2 model (0<a<1, 0<b<3)
param_expow_low <- c(0, 0, 0); param_expow_up <- c(1, 5, 5);
param_hyp1_low <- c(0); param_hyp1_up <- c(1);
param_hyp2_low <- c(0, 0); param_hyp2_up <- c(1, 1);

##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
mle_model_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
mle_model_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp2 <- optim(param2_init, mle_hyp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)

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
p_prd_pow1 <- (t_int+1)^(-parm_pow1[1])
p_prd_pow2 <- parm_pow2[1]*(t_int+1)^(-parm_pow2[2])
p_prd_exp1 <- exp(-parm_exp1[1]*t_int)
p_prd_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*t_int)
p_prd_expow <- parm_expow[1]*exp(-parm_expow[2]*t_int)*(t_int+1)^(-parm_expow[3])
p_prd_hyp1 <- 1/(1+parm_hyp1[1]*t_int)
p_prd_hyp2 <- parm_hyp2[1]/(1+parm_hyp2[2]*t_int)

# Proportion of the explained variances for each model
r2_pow1 = 1-sum((p_corr-p_prd_pow1)^2)/sum((p_corr-mean(p_corr))^2)
r2_pow2 = 1-sum((p_corr-p_prd_pow2)^2)/sum((p_corr-mean(p_corr))^2)
r2_exp1 = 1-sum((p_corr-p_prd_exp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_exp2 = 1-sum((p_corr-p_prd_exp2)^2)/sum((p_corr-mean(p_corr))^2)
r2_expow = 1-sum((p_corr-p_prd_expow)^2)/sum((p_corr-mean(p_corr))^2)
r2_hyp1 = 1-sum((p_corr-p_prd_hyp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_hyp2 = 1-sum((p_corr-p_prd_hyp2)^2)/sum((p_corr-mean(p_corr))^2)

# Generate summary
minus_loglik_MLE = round(c(mle_model_pow1$value, mle_model_pow2$value,
                           mle_model_exp1$value, mle_model_exp2$value, mle_model_expow$value,
                           mle_model_hyp1$value, mle_model_hyp2$value), 3)
r2_mle <- round(c(r2_pow1, r2_pow2, r2_exp1, r2_exp2, r2_expow, r2_hyp1, r2_hyp2), 3)  # round off the value, up to 3 digits
names = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")
pars_mle <- matrix(NA, nrow=3, ncol=7)
pars_mle[,1] <- c(mle_model_pow1$par, NA, NA)
pars_mle[,2] <- c(mle_model_pow2$par, NA)
pars_mle[,3] <- c(mle_model_exp1$par, NA, NA)
pars_mle[,4] <- c(mle_model_exp2$par, NA)
pars_mle[,5] <- mle_model_expow$par
pars_mle[,6] <- c(mle_model_hyp1$par, NA, NA)
pars_mle[,7] <- c(mle_model_hyp2$par, NA)
pars_mle <- round(pars_mle, 3)
dimnames(pars_mle) = list(c('par1', 'par2', 'par3'),c('POW1', 'POW2', 'EXP1', 'EXP2', 'EXPOW', 'HYP1', 'HYP2'))

mle_summary = data.frame(Models = names, loglik = - minus_loglik_MLE, r2 = r2_mle)

##########################
## LSE                  ##
##########################

lse_model_pow1  <- optim(param1_init, lse_pow1,  method="L-BFGS-B", lower=param_pow1_low,  upper=param_pow1_up,  int=t_int, y=p_corr)
lse_model_pow2  <- optim(param2_init, lse_pow2,  method="L-BFGS-B", lower=param_pow2_low,  upper=param_pow2_up,  int=t_int, y=p_corr)
lse_model_exp1  <- optim(param1_init, lse_exp1,  method="L-BFGS-B", lower=param_exp1_low,  upper=param_exp1_up,  int=t_int, y=p_corr)
lse_model_exp2  <- optim(param2_init, lse_exp2,  method="L-BFGS-B", lower=param_exp2_low,  upper=param_exp2_up,  int=t_int, y=p_corr)
lse_model_expow <- optim(param3_init, lse_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, y=p_corr)
lse_model_hyp1  <- optim(param1_init, lse_hyp1,  method="L-BFGS-B", lower=param_hyp1_low,  upper=param_hyp1_up,  int=t_int, y=p_corr)
lse_model_hyp2  <- optim(param2_init, lse_hyp2,  method="L-BFGS-B", lower=param_hyp2_low,  upper=param_hyp2_up,  int=t_int, y=p_corr)

for (i in 1:100) {
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  temp_pow1  <- optim(param1_init, lse_pow1,  method="L-BFGS-B", lower=param_pow1_low,  upper=param_pow1_up,  int=t_int, y=p_corr)
  temp_pow2  <- optim(param2_init, lse_pow2,  method="L-BFGS-B", lower=param_pow2_low,  upper=param_pow2_up,  int=t_int, y=p_corr)
  temp_exp1  <- optim(param1_init, lse_exp1,  method="L-BFGS-B", lower=param_exp1_low,  upper=param_exp1_up,  int=t_int, y=p_corr)
  temp_exp2  <- optim(param2_init, lse_exp2,  method="L-BFGS-B", lower=param_exp2_low,  upper=param_exp2_up,  int=t_int, y=p_corr)
  temp_expow <- optim(param3_init, lse_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, y=p_corr)
  temp_hyp1  <- optim(param1_init, lse_hyp1,  method="L-BFGS-B", lower=param_hyp1_low,  upper=param_hyp1_up,  int=t_int, y=p_corr)
  temp_hyp2  <- optim(param2_init, lse_hyp2,  method="L-BFGS-B", lower=param_hyp2_low,  upper=param_hyp2_up,  int=t_int, y=p_corr)
  
  if(temp_pow1$value < lse_model_pow1$value) lse_model_pow1 <- temp_pow1
  if(temp_pow2$value < lse_model_pow2$value) lse_model_pow2 <- temp_pow2
  if(temp_exp1$value < lse_model_exp1$value) lse_model_exp1 <- temp_exp1
  if(temp_exp2$value < lse_model_exp2$value) lse_model_exp2 <- temp_exp2
  if(temp_expow$value < lse_model_expow$value) lse_model_expow <- temp_expow
  if(temp_hyp1$value < lse_model_hyp1$value) lse_model_hyp1 <- temp_hyp1
  if(temp_hyp2$value < lse_model_hyp2$value) lse_model_hyp2 <- temp_hyp2
}

parm_lse_pow1 <- lse_model_pow1$par
parm_lse_pow2 <- lse_model_pow2$par
parm_lse_exp1 <- lse_model_exp1$par
parm_lse_exp2 <- lse_model_exp2$par
parm_lse_expow <- lse_model_expow$par
parm_lse_hyp1 <- lse_model_hyp1$par
parm_lse_hyp2 <- lse_model_hyp2$par

p_lse_pow1  <- (t_int+1)^(-parm_lse_pow1[1])
p_lse_pow2  <- parm_lse_pow2[1]*(t_int+1)^(-parm_lse_pow2[2])
p_lse_exp1  <- exp(-parm_lse_exp1[1]*t_int)
p_lse_exp2  <- parm_lse_exp2[1]*exp(-parm_lse_exp2[2]*t_int)
p_lse_expow <- parm_lse_expow[1]*exp(-parm_lse_expow[2]*t_int)*(t_int+1)^(-parm_lse_expow[3])
p_lse_hyp1  <- 1/(1+parm_lse_hyp1[1]*t_int)
p_lse_hyp2  <- parm_lse_hyp2[1]/(1+parm_lse_hyp2[2]*t_int)

r2_lse_pow1  <- 1-sum((p_corr-p_lse_pow1)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_pow2  <- 1-sum((p_corr-p_lse_pow2)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_exp1  <- 1-sum((p_corr-p_lse_exp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_exp2  <- 1-sum((p_corr-p_lse_exp2)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_expow <- 1-sum((p_corr-p_lse_expow)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_hyp1  <- 1-sum((p_corr-p_lse_hyp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_lse_hyp2  <- 1-sum((p_corr-p_lse_hyp2)^2)/sum((p_corr-mean(p_corr))^2)

SSE_LSE = round(c(lse_model_pow1$value, lse_model_pow2$value,
                  lse_model_exp1$value, lse_model_exp2$value, lse_model_expow$value,
                  lse_model_hyp1$value, lse_model_hyp2$value), 3)
r2_lse <- round(c(r2_lse_pow1, r2_lse_pow2, r2_lse_exp1, r2_lse_exp2, r2_lse_expow, r2_lse_hyp1, r2_lse_hyp2), 3)  # round off the value, up to 3 digits
pars_lse <- matrix(NA, nrow=3, ncol=7)
pars_lse[,1] <- c(lse_model_pow1$par, NA, NA)
pars_lse[,2] <- c(lse_model_pow2$par, NA)
pars_lse[,3] <- c(lse_model_exp1$par, NA, NA)
pars_lse[,4] <- c(lse_model_exp2$par, NA)
pars_lse[,5] <- lse_model_expow$par
pars_lse[,6] <- c(lse_model_hyp1$par, NA, NA)
pars_lse[,7] <- c(lse_model_hyp2$par, NA)
pars_lse <- round(pars_lse, 3)
dimnames(pars_lse) = list(c('par1', 'par2', 'par3'),c('POW1', 'POW2', 'EXP1', 'EXP2', 'EXPOW', 'HYP1', 'HYP2'))

lse_summary = data.frame(Models = names, SSE = SSE_LSE, r2 = r2_lse)

# Plot the MLE results using simple R graphics. Use ggplots for fancy graphs
x <- seq(0,20, 0.05)

p_pow1 <- (x+1)^(-parm_pow1[1])
p_pow2 <- parm_pow2[1]*(x+1)^(-parm_pow2[2])
p_exp1 <- exp(-parm_exp1[1]*x)
p_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*x)
p_expow <- parm_expow[1]*exp(-parm_expow[2]*x)*(x+1)^(-parm_expow[3])
p_hyp1 <- 1/(1+parm_hyp1[1]*x)
p_hyp2 <- parm_hyp2[1]/(1+parm_hyp2[2]*x)

plot(x, p_pow1, ylim=c(0,1), xlab='Time t', ylab='Proportion Correct', 
     main="MLE results", type='l', lwd=2, col='#e41a1c')
lines(x, p_pow2, lwd=2, lty='dashed', col='#377eb8')
lines(x, p_exp1, lwd=2, col='#4daf4a')
lines(x, p_exp2, lwd=2, col='#984ea3')
lines(x, p_expow, lwd=2, lty='dashed', col='#ff7f00')
lines(x, p_hyp1, lwd=2, col='#a65628')
lines(x, p_hyp2, lwd=2, col='#f781bf')

points(t_int, p_corr, pch=19, cex=1.5)
ltys = c('solid','dashed','solid','solid','dashed','solid','solid')
cols = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf')
legend(14, 1, legend=names, lwd=2, cex=1, lty=ltys, col=cols)

# =========================
# LSE Plot (separate figure)
# =========================

x <- seq(0,20, 0.05)

p_pow1_lse <- (x+1)^(-parm_lse_pow1[1])
p_pow2_lse <- parm_lse_pow2[1]*(x+1)^(-parm_lse_pow2[2])
p_exp1_lse <- exp(-parm_lse_exp1[1]*x)
p_exp2_lse <- parm_lse_exp2[1]*exp(-parm_lse_exp2[2]*x)
p_expow_lse <- parm_lse_expow[1]*exp(-parm_lse_expow[2]*x)*(x+1)^(-parm_lse_expow[3])
p_hyp1_lse <- 1/(1+parm_lse_hyp1[1]*x)
p_hyp2_lse <- parm_lse_hyp2[1]/(1+parm_lse_hyp2[2]*x)

plot(x, p_pow1_lse, ylim=c(0,1), xlab='Time t', ylab='Proportion Correct', 
     main="LSE results", type='l', lwd=2, col='#e41a1c')

lines(x, p_pow2_lse, lwd=2, lty='dashed', col='#377eb8')
lines(x, p_exp1_lse, lwd=2, col='#4daf4a')
lines(x, p_exp2_lse, lwd=2, col='#984ea3')
lines(x, p_expow_lse, lwd=2, lty='dashed', col='#ff7f00')
lines(x, p_hyp1_lse, lwd=2, col='#a65628')
lines(x, p_hyp2_lse, lwd=2, col='#f781bf')

points(t_int, p_corr, pch=19, cex=1.5)

legend(14, 1, legend=names, lwd=2, cex=1, lty=ltys, col=cols)

# print maximized likehood values
print('- MLE results ------------')
print(mle_summary,4)

# print bet-fit parameter values
print('- Best-fit MLE parameters --------')
print(pars_mle,4)

print('- LSE results ------------')
print(lse_summary,4)

print('- Best-fit LSE parameters --------')
print(pars_lse,4)
