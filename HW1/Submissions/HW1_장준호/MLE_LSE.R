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


# POW 1
mle_pow1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- (int + 1)^(-param[1])
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) ) 
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

# POW 2
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

#EXP1
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


# EXP2
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

#EXPOW
mle_expow <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp((-param[2])*int)*(int + 1)^(-param[3])
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

#HYP1
mle_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- 1 / (1 + (param[1] * int))
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

#HYP2
mle_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1] / (1 + (param[2] * int))
  
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1)*( x*log(p) + (n-x)*log(1-p) )
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}




#====================================================#
print('least squares estimates')

# POW 1
mle_pow1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- (int + 1)^(-param[1])
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  

}

# POW 2
mle_pow2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*(1+int)^(-param[2])
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}

#EXP1
mle_exp1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- exp((-param[1])*int)
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}


# EXP2
mle_exp2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp((-param[2])*int)
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}

#EXPOW
mle_expow <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1]*exp((-param[2])*int)*(int + 1)^(-param[3])
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}

#HYP1
mle_hyp1 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- 1 / (1 + (param[1] * int))
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}

#HYP2
mle_hyp2 <- function(param, int=t_int, n=n_total, x=n_corr)  {
  # predicted probability by parameters
  p <- param[1] / (1 + (param[2] * int))
  
  # calculate correct rate
  corr_p <- x/n
  
  # Calculate SSE
  sse <- sum((corr_p - p)^2)
  
  # return sse
  return(sse)  
}
