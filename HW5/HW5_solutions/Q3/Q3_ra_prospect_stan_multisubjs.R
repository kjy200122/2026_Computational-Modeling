# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)

rm(list=ls())  # remove all variables 

library(rstan)

# set seed
set.seed(08826)

# source HDIofMCMC.R to calculate HDI
source("Q3/HDIofMCMC.R") 

# read the data file
dat = read.table("ra_exampleData.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject (=140)
numPars = 3               # number of parameters

# make a list of data
gain     = array(NA, c(N, T))
loss     = array(NA, c(N, T))
cert     = array(NA, c(N, T))
gamble   = array(NA, c(N, T))

for (i in 1:N){
  ID = allSubjs[i]
  tmpData = subset(dat, subjID==ID)
  
  gain[i,] = tmpData$gain
  loss[i,] = abs(tmpData$loss)
  cert[i,] = tmpData$cert
  gamble[i,] = tmpData$gamble
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
output_Q3 = stan("Q3/ra_prospect_multisubjs.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 4000, warmup=1000, chains=4, cores=4)

save(output_Q3, file = 'Q3/output_Q3_ra.RData')

# traceplot
traceplot(output_Q3)

# print summary
print(output_Q3)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output_Q3)

# plot posteriors 
hist(parameters$rho)
hist(parameters$lambda)
hist(parameters$tau)

# 95% HDI of rho
HDIofMCMC(parameters$rho, credMass = 0.95)
HDIofMCMC(parameters$lambda, credMass = 0.95)
HDIofMCMC(parameters$tau, credMass = 0.95)
