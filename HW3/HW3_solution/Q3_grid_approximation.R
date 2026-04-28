####################################
########## Basic Settings ##########
####################################
rm(list=ls())  # clear workspace
graphics.off() # close all figures

set.seed(2024)  # set a seed number for replication

library(tidyverse) 

# Import data
data = read.table(file = './ra_exampleData.txt', sep = '\t', header = TRUE)

subjectList = unique(data$subjID)
N = length(subjectList)
trials = 140

####################################
############# Functions ############
####################################

# Filter data with subjID


###  Function for calculating likelihood of RA prospect model ###
likelihood_ra_prospect <- function(param, T = trials, data_subject) {
  # Prepare data
  gain = data_subject$gain
  loss = data_subject$loss
  cert = data_subject$cert
  gamble = data_subject$gamble
  
  
  # Parameters
  rho = param[1]
  lambda = param[2]
  tau = param[3]

  # Initialize sum of minus log likelihood
  sum_minusLL = 0  
  
  # Calculate minus log likelihood
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL      
  }
  
  # Return the likelihood
  likelihood = exp(-sum_minusLL)
  return(likelihood)
}


### Function for estimating posterior using grid approximation ###
# step 1. Define the grid 
# step 2. compute prior at each parameter value 
# step 3. compute likelihood at each parameter value
# step 4. compute unstandardized posterior at each parameter value
# step 5. Standardize the posterior by dividing by the sum of all values. 

grid_approximation <- function(data_subject, T = trials, n_grid){
  
  # step 1. Define the grid 
  rho_grids = seq(from = 0, to = 2, length.out = n_grid) 
  lambda_grids = seq(from = 0, to = 10, length.out = n_grid) 
  tau_grids = seq(from = 0, to = 5, length.out = n_grid) 
  
  # step 2. compute prior at each parameter value (uniform prior)
  df_prior <- data.frame()
  i = 1
  for (rho in rho_grids){
    for (lambda in lambda_grids){
      for (tau in tau_grids){
        df_prior[i,'rho']    = rho
        df_prior[i,'lambda'] = lambda
        df_prior[i,'tau']    = tau
        df_prior[i, 'prior']  = 1 / n_grid^3
        i  = i + 1
        
      }
    }
  }
  
  
  # step 3. compute likelihood at each parameter value
  df_likelihood <- data.frame()
  i = 1
  for (rho in rho_grids){
    for (lambda in lambda_grids){
      for (tau in tau_grids){
        df_likelihood[i,'rho']    = rho
        df_likelihood[i,'lambda'] = lambda
        df_likelihood[i,'tau']    = tau
        df_likelihood[i, 'likelihood']  = likelihood_ra_prospect(c(rho, lambda, tau), T = trials, data_subject) 
        i  = i + 1
      }
    }
  }    
  
  # step 4. compute unstandardized posterior at each parameter value
  df_posterior <- data.frame()
  i = 1
  for (rho in rho_grids){
    for (lambda in lambda_grids){
      for (tau in tau_grids){
        df_posterior[i,'rho']    = rho
        df_posterior[i,'lambda'] = lambda
        df_posterior[i,'tau']    = tau
        df_posterior[i, 'posterior']  = df_prior[i,'prior'] * df_likelihood[i,'likelihood']
        i  = i + 1
      }
    }
  }    

  # step 5. Standardize the posterior by dividing by the sum of all values.
  df_posterior$posterior = df_posterior$posterior / sum(df_posterior$posterior)
  
  return(df_posterior)
}


# grid approximation with 10 grids
n_grid = 10
posteriors_10_grids <- list()
for (s in 1:length(subjectList)){
  # subject ID
  subject <- subjectList[s]
  
  # Prepare data
  data_subject = data[data$subjID == subject,]

  # Run 
  posteriors_10_grids[[s]] <- grid_approximation(data_subject, T = trials, n_grid)
}

# Save data
save(posteriors_10_grids, file = "posteriors_10_grids.RData")


# grid approximation with 30 grids
n_grid = 30 
posteriors_30_grids <- list()
for (s in 1:length(subjectList)){
  # subject ID
  subject <- subjectList[s]
  
  # Prepare data
  data_subject = data[data$subjID == subject,]
  
  # Run 
  posteriors_30_grids[[s]] <- grid_approximation(data_subject, T = trials, n_grid)
}

# Save data
save(posteriors_30_grids, file = "posteriors_30_grids.RData")

  


