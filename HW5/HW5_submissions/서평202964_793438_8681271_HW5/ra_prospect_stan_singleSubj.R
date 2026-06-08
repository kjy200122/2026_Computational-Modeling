# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)

rm(list=ls())  # remove all variables 

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

dat_all = dat[, c("subjID", "gain", "loss", "cert", "gamble")]

#Make a matrix for each data
gain_mat   <- matrix(NA, nrow=N, ncol=T)
loss_mat   <- matrix(NA, nrow=N, ncol=T)
cert_mat   <- matrix(NA, nrow=N, ncol=T)
gamble_mat <- matrix(NA, nrow=N, ncol=T)


# use first subject only
for (n in 1:N){

    tmpData = subset(dat, subjID==allSubjs[n])
    gain_mat[n, ]   <- tmpData$gain
    loss_mat[n, ]   <- abs(tmpData$loss) #absolute value
    cert_mat[n, ]   <- tmpData$cert
    gamble_mat[n, ] <- tmpData$gamble
    
}
#save to datalist
dataList <- list(
  N = N,
  T = T,
  gain = gain_mat,
  loss = loss_mat,
  cert = cert_mat,
  gamble = gamble_mat)

# run!
output = stan("ra_prospect_singleSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 4000, warmup=1000, chains=4, cores=4)

# traceplot
rstan::traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

# plot posteriors 
hist(parameters$rho)
hist(parameters$lambda)
hist(parameters$tau)

# 95% HDI of rho
HDIofMCMC(parameters$lambda, credMass = 0.95)

#Means of posteriors for each parameter
rho_mean <- apply(parameters$rho, 2, mean)
lambda_mean <- apply(parameters$lambda, 2, mean)
tau_mean <- apply(parameters$tau, 2, mean)

posterior_means <- data.frame(
  subjID = allSubjs,
  rho_mean = rho_mean,
  lambda_mean = lambda_mean,
  tau_mean = tau_mean
)

posterior_means
