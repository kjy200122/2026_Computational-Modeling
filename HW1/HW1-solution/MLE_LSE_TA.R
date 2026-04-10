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

################################################################################
# Functions for MLE & LSE
################################################################################

mle_function <- function(func){
  ret <- function(param, int=t_int, n=n_total, x=n_corr)  {
    # predicted probability by parameters
    p <- func(int,param) 
    
    # ensure 0 < p < 1
    p[p<=0] <- 10e-6
    p[p>=1] <- 1 - 10e-6
    
    # Calculate minus log-likelihood
    loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) ) 
    
    # return the summed minus log-likelihood as the function value
    sum(loglik)
  }
  return(ret)
}

lse_function <- function(func){
  ret <- function(param, int=t_int, n=n_total, x=n_corr)  {
    # predicted probability by parameters
    p <- func(int,param) 
    
    # ensure 0 < p < 1
    p[p<=0] <- 10e-6
    p[p>=1] <- 1 - 10e-6
    
    p_obs <- x/n
    
    # Calculate the sum of squared errors
    sse <- sum((p - p_obs) ** 2)    
    return(sse)
  }
  return(ret)
}

################################################################################
# Functions to compute the probability
################################################################################

compute_p_pow1  <- function(int, param) { (1 + int) ** (-param[1]) }
compute_p_pow2  <- function(int, param) { param[1] * (1 + int) ** (-param[2]) }
compute_p_exp1  <- function(int, param) { exp((-param[1]) * int) }
compute_p_exp2  <- function(int, param) { param[1] * exp((-param[2]) * int) }
compute_p_expow <- function(int, param) { param[1] * exp((-param[2]) * int) * (1 + int) ** (-param[3]) }
compute_p_hyp1  <- function(int, param) { 1 / (1 + param[1] * int) }
compute_p_hyp2  <- function(int, param) { param[1] / (1 + param[2] * int) }

################################################################################
# Loss functions for MLE
################################################################################

mle_pow1  <- mle_function(compute_p_pow1)
mle_pow2  <- mle_function(compute_p_pow2)
mle_exp1  <- mle_function(compute_p_exp1)
mle_exp2  <- mle_function(compute_p_exp2)
mle_expow <- mle_function(compute_p_expow)
mle_hyp1  <- mle_function(compute_p_hyp1)
mle_hyp2  <- mle_function(compute_p_hyp2)

################################################################################
# Loss functions for LSE
################################################################################

lse_pow1  <- lse_function(compute_p_pow1)
lse_pow2  <- lse_function(compute_p_pow2)
lse_exp1  <- lse_function(compute_p_exp1)
lse_exp2  <- lse_function(compute_p_exp2)
lse_expow <- lse_function(compute_p_expow)
lse_hyp1  <- lse_function(compute_p_hyp1)
lse_hyp2  <- lse_function(compute_p_hyp2)
