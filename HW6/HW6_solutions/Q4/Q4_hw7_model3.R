# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
# The only difference between model1 is data & stan file name

#rm(list=ls())  # remove all variables

library(rstan)

# read the data file
dat = read.table("simul_data_hw7_model3.txt", header=T, sep="\t")

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
output_Q4 = stan("Q4/Q4_hw7_model3.stan", 
                 data = dataList, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta", "sigma", 
                                           "alpha_pos", "alpha_neg", "beta"),
                 iter = 2000, warmup=1000, chains=2, cores=2)

save(output_Q4, file = "Q4/output_Q4_model3.RData")
