rm(list=ls())

ra_LL_given_param <- function(param_grid, n_trials, data){
  # param grid is list consisted with vector: length(param_grid$parameter) = # of grid
  # Thus, output log_lik is also a vector - length: (# of grid)
  log_lik = rep(0, nrow(param_grid))
  rho = param_grid$rho
  lambda = param_grid$lambda
  tau = param_grid$tau
  
  gain = data$gain
  loss = data$loss
  cert = data$cert
  gamble = data$gamble
  
  for (t in 1:n_trials) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_log_lik = log(pGamble)*gamble[t] + log(1-pGamble)*(1-gamble[t])
    log_lik = log_lik + tmp_log_lik
  }
  return(log_lik)
}

ra_grid_approx <- function(data, len_grid, param_grid, log_prior){
  # len_grid: the number of gird "per parameter"
  # param_grid: the entire combinations of parameter candidates
  # log prior: list consisted with functions specifying the prior distribution for each parameter
  
  n_sub = length(unique(data$subjID))
  
  rho_all_sub = matrix(data=0, nrow=n_sub, ncol=len_grid)
  lambda_all_sub = matrix(data=0, nrow=n_sub, ncol=len_grid)
  tau_all_sub= matrix(data=0, nrow=n_sub, ncol=len_grid)
  
  s=1
  for (subjID in unique(data$subjID)){
    
    data_sub = data[data$subjID == subjID, ]
    n_trials = length(data_sub$gamble)
    
    # prior of a certain grid point[i,j,k] = p(rho_i)*p(lambda_j)*p(tau_k)
    log_prior_grid = (log_prior$rho(param_grid$rho) 
                      + log_prior$lambda(param_grid$lambda)
                      + log_prior$tau(param_grid$tau))
    
    # get log liklihood for each grid point
    log_lik_grid = ra_LL_given_param(param_grid, n_trials, data_sub)
    
    # get p(theta[i,j,k] | data) * p(theta[i,j,k]) (proportionate to p(theta[i,j,k]))
    log_pos_grid = log_lik_grid+log_prior_grid
    
    rho_dist = rep(0, len_grid)
    lambda_dist = rep(0, len_grid)
    tau_dist = rep(0, len_grid)
    
    for (i in 1:len_grid){
      
      # Fixed one parameter, sum across the all grid spaces
      rho_dist[i] = sum(exp(log_pos_grid[param_grid$rho == rho_grid[i]]))
      lambda_dist[i] = sum(exp(log_pos_grid[param_grid$lambda == lambda_grid[i]]))
      tau_dist[i] = sum(exp(log_pos_grid[param_grid$tau == tau_grid[i]]))
    }
    
    rho_all_sub[s, ] = rho_dist/sum(rho_dist)
    lambda_all_sub[s, ]= lambda_dist/sum(lambda_dist)
    tau_all_sub[s, ] = tau_dist/sum(tau_dist)
    s = s+1
  }
  return(list(rho=rho_all_sub, lambda=lambda_all_sub, tau=tau_all_sub))
}


data=read.delim("ra_exampleData.txt")

# prior: assume uniform prior (with log)
# another prior can be tested, by revising this function
rho_prior <- function(x){
  return(dunif(x=x, min=0, max=2, log=TRUE))
}

lambda_prior <- function(x){
  return(dunif(x=x, min=0, max=10, log=TRUE))
}

tau_prior <- function(x){
  return(dunif(x=x, min=0, max=5, log=TRUE))
}

### Test 10 grids per parameter
len_grid = 10
rho_grid = seq(0.01, 1.99, length.out = len_grid)
lambda_grid = seq(0.01, 4.99, length.out = len_grid)
tau_grid = seq(0.01, 9.99, length.out = len_grid)
param_grid = expand.grid(rho=rho_grid, lambda=lambda_grid, tau=tau_grid)

log_prior = list(rho=rho_prior, lambda=lambda_prior, tau=tau_prior)
result = ra_grid_approx(data=data, len_grid=len_grid, param_grid = param_grid, log_prior=log_prior)

par(mfrow=c(5,3))

for (s in 1:nrow(result$rho)){
  
 plot(
    rho_grid,
    result$rho[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="rho",
    ylab="p(rho | data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
  
 plot(
    lambda_grid,
    result$lambda[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="lambda",
    ylab="p(lambda | data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
  
  plot(
    tau_grid,
    result$tau[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="tau",
    ylab="p(tau | data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
}

## Test 30 grids per parameter
len_grid = 30
rho_grid = seq(0.01, 1.99, length.out = len_grid)
lambda_grid = seq(0.01, 4.99, length.out = len_grid)
tau_grid = seq(0.01, 9.99, length.out = len_grid)
param_grid = expand.grid(rho=rho_grid, lambda=lambda_grid, tau=tau_grid)

log_prior = list(rho=rho_prior, lambda=lambda_prior, tau=tau_prior)
result = ra_grid_approx(data=data, len_grid=len_grid, param_grid = param_grid, log_prior=log_prior)

par(mfrow=c(5,3))

for (s in 1:nrow(result$rho)){
  
  plot(
    rho_grid,
    result$rho[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="rho",
    ylab="p(rho|data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
  
  plot(
    lambda_grid,
    result$lambda[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="lambda",
    ylab="p(lambda|data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
  
  plot(
    tau_grid,
    result$tau[s,],
    type="l",
    main = paste("ID:", unique(data$subjID)[s]),
    xlab="tau",
    ylab="p(tau|data)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
  )
}

