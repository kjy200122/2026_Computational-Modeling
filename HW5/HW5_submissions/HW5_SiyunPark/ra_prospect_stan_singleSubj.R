# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)

rm(list=ls())  # remove all variables 

library(rstan)

# set seed
set.seed(08826)

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 
 
# read the data file
data = read.table("ra_exampleData.txt", header=T, sep="\t")
subjs <- c(2, 3, 4, 6, 7)

# allSubjs = unique(dat$subjID)  # all subject IDs
N = length(subjs)      # number of subjects
T = table(data$subjID)[1]  # number of trials per subject (=140)
#numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

gain <- matrix(NA, nrow = N, ncol = T)
loss <- matrix(NA, nrow = N, ncol = T)
cert <- matrix(NA, nrow = N, ncol = T)
gamble <- matrix(NA, nrow = N, ncol = T)

for (s in 1:N) {
  tmpData <- subset(data, subjID == subjs[s])
  
  gain[s, ] <- tmpData$gain
  loss[s, ] <- abs(tmpData$loss)
  cert[s, ] <- tmpData$cert
  gamble[s, ] <- tmpData$gamble
}

dataList <- list(
  N = N,
  T = T,
  gain = gain,
  loss = loss,
  cert = cert,
  gamble = gamble
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

# plot posteriors 
par(mfrow = c(5, 3))

for (s in 1:N) {
  hist(parameters$rho[, s], main = paste("Subject", subjs[s], "rho"), xlab = "rho")
  hist(parameters$lambda[, s], main = paste("Subject", subjs[s], "lambda"), xlab = "lambda")
  hist(parameters$tau[, s], main = paste("Subject", subjs[s], "tau"), xlab = "tau")
}

# 95% HDI of rho
# HDIofMCMC(parameters$lambda, credMass = 0.95)

posterior_means <- data.frame(
  subjID = subjs,
  rho = apply(parameters$rho, 2, mean),
  lambda = apply(parameters$lambda, 2, mean),
  tau = apply(parameters$tau, 2, mean)
)

print(posterior_means)
