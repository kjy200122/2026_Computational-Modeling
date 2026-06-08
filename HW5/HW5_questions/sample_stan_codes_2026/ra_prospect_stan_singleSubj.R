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

# use first subject only
tmpData = subset(dat, subjID==2)
dataList <- list(
  T       = T,
  gain    = tmpData$gain,
  loss    = abs(tmpData$loss),   # absolute value
  cert    = tmpData$cert,
  gamble  = tmpData$gamble
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
hist(parameters$rho)
hist(parameters$lambda)
hist(parameters$tau)

# 95% HDI of rho
HDIofMCMC(parameters$lambda, credMass = 0.95)
