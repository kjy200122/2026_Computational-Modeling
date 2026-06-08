# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
# The only difference between model1 is data & stan file name

#rm(list=ls())  # remove all variables

library(rstan)

# read the data file
dat = read.table("simul_data_hw7_model2.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject 

choice  <- array(-1, c(N, T))
outcome <- array(0, c(N, T))

for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T] <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}

dataList <- list(
  N       = N,
  T       = T,
  Tsubj   = rep(T, N),
  choice  = choice,
  outcome = outcome
)

# run!
output_Q3 = stan("Q3/Q3_hw7_model2.stan", data = dataList, pars = c("mu_alpha", "mu_beta", "sigma", "alpha", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

save(output_Q3, file = "Q3/output_Q3.RData")
