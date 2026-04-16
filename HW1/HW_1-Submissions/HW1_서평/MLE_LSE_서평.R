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
mle_pow1 <-  function(param, int=t_int, n=n_total, x=n_corr)  {
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

mle_exp1 <-function(param, int=t_int, n=n_total, x=n_corr)  {
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
  p <- param[1]*exp((-param[2])*int)*((int+1)^(-param[3]))
  
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
  p <- 1/(1 + (param[1]*int))
  
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
  p <- param[1]/(1 + (param[2]*int))
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}


#LSE (Least Squares Estimation) Functions 

lse_pow1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n  # 실제 관측 비율 계산
  p_pred <- (1+int)^(-param[1])
  
  # p_pred를 0~1 범위로 제한, 무한대로 가는 오류를 막는다
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  # LSE 계산 : 잔차의 제곱의 합
  sum((p_obs - p_pred)^2)
}

lse_pow2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- param[1]*(1+int)^(-param[2])
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}

lse_exp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- exp((-param[1])*int)
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}

lse_exp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- param[1]*exp((-param[2])*int)
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}

lse_expow <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- param[1]*exp((-param[2])*int)*((int+1)^(-param[3]))
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}

lse_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- 1/(1 + (param[1]*int))
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}

lse_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr) {
  p_obs <- x / n
  p_pred <- param[1]/(1 + (param[2]*int))
  
  p_pred[p_pred <= 0] <- 1e-6
  p_pred[p_pred >= 1] <- 1 - 1e-6
  
  sum((p_obs - p_pred)^2)
}
