
# library hBayesDM to get ra_prospect model
library(hBayesDM)

# read ra_exampleData
ra_example <- read.table("ra_exampleData.txt", header=TRUE)
subj_list <- c(2, 3, 4, 6, 7)

output <- ra_prospect(
  data = "ra_exampleData.txt",
  niter = 4000,
  nwarmup = 1000,
  nchain = 4,
  ncore = 2,
  nthin = 1,
  inits = "vb",
  indPars = "mean",
  modelRegressor = FALSE,
  vb = FALSE,
  inc_postpred = FALSE,
  adapt_delta = 0.95,
  stepsize = 1,
  max_treedepth = 10
)

# Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
plot(output, type = "trace")

# Check Rhat values (all Rhat values should be less than or equal to 1.1)
rhat(output)

# Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
plot(output)

# Show the WAIC and LOOIC model fit estimates
printFit(output)

#plot individual posterior distribution
plotInd(output, "rho")
plotInd(output, "lambda")
plotInd(output, "tau")
plotInd(output, c("rho", "lambda", "tau"), show_density = F)

post_mean <- output$allIndPars
post_mean

#To compare the parameter values, do the MLE

N = 5  # number of subjects
T = 140 # number of trials per subject

# cert[t]: outcome of a certain option on trial t
# gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
# gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.
# 
# evSafe: expected value of a certain (safe) option
# evGamble: expected value of a risky option (gamble)
# pGamble   # probability of choosing a gamble on each trial
# free parameters: rho, tau, lambda

ra_prospect_MLE <- function(param, gain, loss, cert, gamble){
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  sum_minusLL
}

param_ra_prospect_low <- c(0, 0, 0); param_ra_prospect_up <- c(2, 10, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)

##########################
## MLE                  ##
##########################

# for all subjects
parm_ra_prospect_all <- data.frame(subjID=subj_list, rho=NA, lambda=NA, tau=NA, minusLL=NA)
for(n in 1:N){
  subj_data <- subset(ra_example, subjID==subj_list[n])
  param3_init <- runif(3)
  
  mle_ra_prospect_all <- optim(param3_init, ra_prospect_MLE, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
  # Try many different inits to escape from the local maxima
  for (i in 1:100)
  {
    # Re-generate random inits. Is it the best way to do this?
    param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
    
    # Do the MLE again
    temp_ra_prospect <- optim(param3_init, ra_prospect_MLE, method="L-BFGS-B", lower=param_ra_prospect_low, upper=param_ra_prospect_up, gain=subj_data$gain, loss=subj_data$loss, cert=subj_data$cert, gamble=subj_data$gamble)
    
    # Replace the results if the latest optimization yields better result

      if(temp_ra_prospect$value < mle_ra_prospect_all$value) mle_ra_prospect_all <- temp_ra_prospect}

    # Save the MLE parameter estimates
  parm_ra_prospect_all$rho[n] <- mle_ra_prospect_all$par[1]
  parm_ra_prospect_all$lambda[n] <- mle_ra_prospect_all$par[2]
  parm_ra_prospect_all$tau[n] <- mle_ra_prospect_all$par[3]
  parm_ra_prospect_all$minusLL[n] <- mle_ra_prospect_all$value
}
print(round(parm_ra_prospect_all, 4)) #all subjects parameters and minusLL


# Plot the MLE and hBayesDM results using simple R graphics

# Merge MLE estimates and posterior means by subject ID
compare_df <- merge(
  parm_ra_prospect_all,
  post_mean,
  by = "subjID",
  suffixes = c("_MLE", "_hBayesDM")
)
print(round(compare_df, 4))

# Plot MLE estimates on x-axis and hBayesDM posterior means on y-axis
plot(compare_df$rho_MLE, compare_df$rho_hBayesDM,
     xlab = "MLE estimate",
     ylab = "Posterior mean",
     main = "rho", pch = 19)
# To compare two values, add x=y line
abline(0, 1, lty = 2)
# Add subject IDs
text(compare_df$rho_MLE, compare_df$rho_hBayesDM,
     labels = compare_df$subjID, pos = 3)

plot(compare_df$lambda_MLE, compare_df$lambda_hBayesDM,
     xlab = "MLE estimate",
     ylab = "Posterior mean",
     main = "lambda", pch = 19, col = 'blue')
# To compare two values, add x=y line
abline(0, 1, lty = 2)
# Add subject IDs under the dots
text(compare_df$lambda_MLE, compare_df$lambda_hBayesDM,
     labels = compare_df$subjID, pos = 3, col = 'blue')

plot(compare_df$tau_MLE, compare_df$tau_hBayesDM,
     xlab = "MLE estimate",
     ylab = "Posterior mean",
     main = "tau", pch = 19, col = 'red')
# To compare two values, add x=y line
abline(0, 1, lty = 2)
# Add subject IDs under the dots
text(compare_df$tau_MLE, compare_df$tau_hBayesDM,
     labels = compare_df$subjID, pos = 3, col = 'red')

