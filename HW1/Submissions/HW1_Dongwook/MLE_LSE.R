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

## MLE Functions ##
mle_pow1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- (int+1)^(-param[1])
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) ) 
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

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

mle_expow <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp((-param[2])*int)*(int+1)^(-param[3])
  
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
  p <- (1+param[1]*int)^(-1)
  
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
  p <- (param[1])*(1+param[2]*int)^(-1)
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}


## LSE Functions ##
# ∑( pᵢ,obs - p(tᵢ) )²
# pᵢ,obs = observed proportions -> p_corr (n_corr/n_total)
# p(tᵢ) = model prediction proportions -> calculated each model

lse_pow1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- (1+int)^(-param[1])
  sum((x/n - p)^2)
}

lse_pow2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- param[1]*(1+int)^(-param[2])
  sum((x/n - p)^2)
}

lse_exp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- exp(-param[1]*int)
  sum((x/n - p)^2)
}

lse_exp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- param[1]*exp(-param[2]*int)
  sum((x/n - p)^2)
}

lse_expow <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- param[1]*exp(-param[2]*int)*(1+int)^(-param[3])
  sum((x/n - p)^2)
}

lse_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- (1+param[1]*int)^(-1)
  sum((x/n - p)^2)
}

lse_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p <- param[1]*(1+param[2]*int)^(-1)
  sum((x/n - p)^2)
}

