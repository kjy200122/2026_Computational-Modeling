#######################################################
## MLE_LSE.R                                         ##
## Maximum Likelihood Estimation (Myung, 2003)       ##
## By Yun Tang, Psychology, OSU                      ##
##                                                   ##
## Functions Calculating Log-Likelihood              ##   
## Code Written on 12/18/2009                        ##
##                                                   ##
## Modified by Joonsuk Park on Jan 28 2015           ##
## Modified by Jay Myung in Feb 2017                 ##
## Modified by Woo-Young Ahn in March 2018           ##
#######################################################
 
mle_pow2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*(1+int)^(-param[2])
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
   
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) ) 
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_exp2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp((-param[2])*int)
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_pow1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- (1+int)^(-param[1])
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_exp1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- exp((-param[1])*int)
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_expow <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp(-param[2]*int)*((int+1)^(-param[3]))
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- 1/(1+param[1]*int)
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

mle_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]/(1+param[2]*int)
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

################################
#functions for LSE 

lse_pow2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: POW2 model
  p_prd <- param[1] * (1 + int)^(-param[2])
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}

lse_exp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: EXP2 model
  p_prd <- param[1] * exp((-param[2]) * int)
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}

lse_pow1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: POW1 model
  p_prd <- (1 + int)^(-param[1])
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}

lse_exp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: EXP1 model
  p_prd <- exp((-param[1]) * int)
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}

lse_expow <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: EXPOW model
  p_prd <- param[1] * exp(-param[2] * int) * ((int + 1)^(-param[3]))
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}

lse_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: HYP1 model
  p_prd <- 1 / (1 + param[1] * int)
  
  sum((p_obs - p_prd)^2)
}

lse_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  # observed proportion correct
  p_obs <- x / n
  
  # predicted probability by parameters: HYP2 model
  p_prd <- param[1] / (1 + param[2] * int)
  
  # return the summed squared error as the function value
  sum((p_obs - p_prd)^2)
}
