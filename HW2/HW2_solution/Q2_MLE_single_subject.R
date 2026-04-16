rm(list=ls())  # clear workspace
graphics.off() # close all figures

set.seed(202604)  # set a seed number for replication
source('ra_MLE_prototype.R')

# Loading data
data = read.table("ra_exampleData.txt", header = TRUE)

N = 5  # number of subjects
trial = 140 # number of trials per subject
subjList = unique(data$subjID)


k_ra_prospect <- 3 # number of parameters

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param_init <- runif(3)
param_low <- c(0, 0, 0); param_up <- c(2, 10, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)


# data of each subjects
subject = 2 # put subject number that you are interested 
subdata = data[data$subjID == subject,]


##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_prospect <- optim(param_init, mle_ra_prospect, method = "L-BFGS-B", lower = param_low, upper = param_up, T = trial, cert = subdata$cert, gain = subdata$gain, loss = subdata$loss, gamble = subdata$gamble)

# Try many different inits to escape from the local maxima
for (i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param_init <- runif(3); 
  
  # Do the MLE again
  temp_ra_prospect  <- optim(param_init, mle_ra_prospect, lower = param_low, upper = param_up, method = 'L-BFGS-B',
                             T = trial, cert = subdata$cert, gain = subdata$gain, loss = subdata$loss, gamble = subdata$gamble)
  
  # Replace the results if the latest optimization yields better result
  if (temp_ra_prospect$value < mle_model_prospect$value) mle_model_prospect <- temp_ra_prospect
}

# Save the MLE parameter estimates
parm_ra_prospect <- mle_model_prospect$par

df_estimates_single <- data.frame(t(parm_ra_prospect))
colnames(df_estimates_single) <- c('rho', 'lambda', 'tau')

                                  
