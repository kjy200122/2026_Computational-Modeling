# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)

rm(list=ls())  # remove all variables 

getwd()
list.files()

library(rstan)

# set seed
set.seed(08826)

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 
 
# read the data file
dat = read.table("ra_exampleData.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject (=140)
#numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

gain <- matrix(NA, nrow = N, ncol = T)
loss <- matrix(NA, nrow = N, ncol = T)
cert <- matrix(NA, nrow = N, ncol = T)
gamble <- matrix(NA, nrow = N, ncol = T)

# Create matrices for input
for (s in 1:N){
  subj <- allSubjs[s]
  dat_s = dat[dat$subjID == subj, ]
  
  for (j in 1:nrow(dat_s)){
  gain[s,j] = dat_s$gain[j]
  loss[s,j] = abs(dat_s$loss[j])
  cert[s,j] = dat_s$cert[j]
  gamble[s,j] = dat_s$gamble[j]
  }
}

range(loss)

# use all subjects
dataList <- list(
  N = N,
  T = T,
  gain    = gain,
  loss    = loss,   # absolute value
  cert    = cert,
  gamble  = gamble
)

# run!
output = stan("ra_prospect_singleSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 2000, warmup=1000, chains=4, cores=4)

# traceplot
rstan::traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

result_df <- data.frame(
  subjID = allSubjs,
  rho_mean = NA,
  lambda_mean = NA,
  tau_mean = NA
)
for (i in 1:N){
  result_df$rho_mean[i] = mean(parameters$rho[, i])
  result_df$lambda_mean[i] = mean(parameters$lambda[, i])
  result_df$tau_mean[i] = mean(parameters$tau[, i])
}
result_df

# plot posteriors 
hist(parameters$rho[, 1])
hist(parameters$lambda[, 1])
hist(parameters$tau[, 1])

# 95% HDI of rho
HDIofMCMC(parameters$lambda[, 1], credMass = 0.95)
