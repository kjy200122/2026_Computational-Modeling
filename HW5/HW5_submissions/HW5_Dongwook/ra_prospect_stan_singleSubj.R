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

gamble_mat <- rbind(
  subset(dat, subjID == 2)$gamble,
  subset(dat, subjID == 3)$gamble,
  subset(dat, subjID == 4)$gamble,
  subset(dat, subjID == 6)$gamble,
  subset(dat, subjID == 7)$gamble
)

cert_mat <- rbind(
  subset(dat, subjID == 2)$cert,
  subset(dat, subjID == 3)$cert,
  subset(dat, subjID == 4)$cert,
  subset(dat, subjID == 6)$cert,
  subset(dat, subjID == 7)$cert
)

gain_mat <- rbind(
  subset(dat, subjID == 2)$gain,
  subset(dat, subjID == 3)$gain,
  subset(dat, subjID == 4)$gain,
  subset(dat, subjID == 6)$gain,
  subset(dat, subjID == 7)$gain
)

loss_mat <- rbind(
  abs(subset(dat, subjID == 2)$loss),
  abs(subset(dat, subjID == 3)$loss),
  abs(subset(dat, subjID == 4)$loss),
  abs(subset(dat, subjID == 6)$loss),
  abs(subset(dat, subjID == 7)$loss)
)

dataList <- list(
  N      = N,
  T      = T,
  gain   = gain_mat,
  loss   = loss_mat,
  cert   = cert_mat,
  gamble = gamble_mat
)

# run!
output = stan("ra_prospect_singleSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 4000, warmup=1000, chains=4, cores=4)

# traceplot
rstan::traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

# plot posteriors means
# subj 2
mean(parameters$rho[, 1])
mean(parameters$lambda[, 1])
mean(parameters$tau[, 1])
# subj 3
mean(parameters$rho[, 2])
mean(parameters$lambda[, 2])
mean(parameters$tau[, 2])
# subj 4
mean(parameters$rho[, 3])
mean(parameters$lambda[, 3])
mean(parameters$tau[, 3])
# subj 6
mean(parameters$rho[, 4])
mean(parameters$lambda[, 4])
mean(parameters$tau[, 4])
# subj 7
mean(parameters$rho[, 5])
mean(parameters$lambda[, 5])
mean(parameters$tau[, 5])

