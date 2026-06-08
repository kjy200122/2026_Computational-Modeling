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

dat = dat[order(dat$subjID),] # assure data ordered by subjID

#numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

#reorganize data as Stan expects: 700 -> (5(N), 140(T))
gain <- matrix(dat$gain, nrow = N, ncol = T, byrow = TRUE) 
loss <- matrix(dat$loss, nrow = N, ncol = T, byrow = TRUE)
cert <- matrix(dat$cert, nrow = N, ncol = T, byrow = TRUE)
gamble <- matrix(dat$gamble, nrow = N, ncol = T, byrow = TRUE)

dataList <- list(
  N       = N,
  T       = T,
  gain    = gain,
  loss    = abs(loss),   # absolute value
  cert    = cert,
  gamble  = gamble
)

# run!
output = stan("ra_prospect_singleSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
              iter = 2000, warmup=1000, chains=4, cores=4)

# traceplot
rstan::traceplot(output, pars="rho")
rstan::traceplot(output, pars="lambda")
rstan::traceplot(output, pars="tau")

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)
# plot posteriors 
par(mfrow = c(N, 3))
for (i in 1:N) {
  s = unique(dat$subjID)[i]
  
  hist(parameters$rho[, i],
       breaks=30,
       main=sprintf("rho[id=%d]", s),
       xlab="rho",
       xlim=c(0.5,1.5))
  
  hist(parameters$lambda[, i],
       breaks=30,
       main=sprintf("lambda[id=%d]", s),
       xlab="lambda",
       xlim=c(0.5, 3.5))
  
  hist(parameters$tau[, i],
       breaks=30,
       main=sprintf("tau[id=%d]", s),
       xlab="tau",
       xlim=c(0, 7))
}

# calculate posterior means and 95% HDI of parameters of 5 subjects using "apply" for each column(subject)
param_dist_df = data.frame(
  subjID = unique(dat$subjID),
  rho_mean = round(apply(parameters$rho, 2, mean),3),
  rho_2_5 = round(apply(parameters$rho, 2, HDIofMCMC)[1, ],3),
  rho_97_5 = round(apply(parameters$rho, 2, HDIofMCMC)[2, ],3),
  lambda_mean = round(apply(parameters$lambda, 2, mean),3),
  lambda_2_5 = round(apply(parameters$lambda, 2, HDIofMCMC)[1, ],3),
  lambda_97_5 = round(apply(parameters$lambda, 2, HDIofMCMC)[2, ],3),
  tau_mean = round(apply(parameters$tau, 2, mean),3),
  tau_2_5 = round(apply(parameters$tau, 2, HDIofMCMC)[1, ],3),
  tau_97_5 = round(apply(parameters$tau, 2, HDIofMCMC)[2, ],3)
)
print(param_dist_df)
write.csv(param_dist_df, "ra_prospect_posterior_dist.csv")
