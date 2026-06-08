# ra_prospect_stan_allSubj.R
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

# use all subjects
gain<- matrix(NA, nrow = N, ncol = T) # make matrix of gain, loss, cert, gamble
loss <- matrix(NA, nrow = N, ncol = T)
cert <- matrix(NA, nrow = N, ncol = T)
gamble <- matrix(NA, nrow = N, ncol = T)

for (n in 1:N) { #use for loop for all subjects
tmpData = subset(dat, subjID==allSubjs[n])
 gain[n, ]    = tmpData$gain
 loss[n, ]    = abs(tmpData$loss)   # absolute value
 cert[n, ]    = tmpData$cert
gamble[n, ]  = tmpData$gamble
}

dataList <- list(
  N       = N,
  T       = T,
  gain = gain,
  loss = loss,   # absolute value
  cert = cert,
  gamble = gamble
)

# run!
output = stan("ra_prospect_allSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 2000, warmup=1000, chains=4, cores=4)

# traceplot of each parameter
rstan::traceplot(output, pars = "rho")
rstan::traceplot(output, pars = "lambda")
rstan::traceplot(output, pars = "tau")

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

# plot posteriors of all subjects
for (n in 1:N) {
hist(parameters$rho[ ,n], main = paste("Subject", allSubjs[n], "rho"))
hist(parameters$lambda[ ,n], main = paste("Subject", allSubjs[n], "lambda"))
hist(parameters$tau[ ,n], main = paste("Subject", allSubjs[n], "tau"))
}

# 95% HDI of rho
for (n in 1:N) {
HDIofMCMC(parameters$lambda[ ,n], credMass = 0.95)
  print(HDIofMCMC(parameters$lambda[ ,n], credMass = 0.95))}

# get posterior mean of subjects
post_mean <- data.frame(
  subjID = allSubjs,
  rho = colMeans(parameters$rho),
  lambda = colMeans(parameters$lambda),
  tau = colMeans(parameters$tau)
)

print(round(post_mean, 4))
