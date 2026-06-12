# Q1_model1_looic.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)

rm(list=ls())  # remove all variables
# install.packages("loo")  # install loo first if not installed

library(rstan)
library(loo)  

# read the data file
dat = read.table("simul_data_model1.txt", header=T, sep="\t")

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
output = stan("model1_looic.stan", data = dataList, pars = c("mu_alpha", "mu_beta", "sigma", "alpha", "beta", 
                                                                   "log_lik"),
              iter = 2000, warmup=1000, chains=4, cores=8, seed = 2026)
save(output, file = "output_looic.RData")

# m = stan_model("hw5_model1.stan")
# output_vb = vb(m, data = dataList, pars = c("alpha", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),
# )


#####################################
load("output_looic.RData")

Q1_log_lik <- loo::extract_log_lik(stanfit = output, 
                                   parameter_name = "log_lik",
                                   merge_chains = FALSE)
Q1_r_eff <- relative_eff(exp(Q1_log_lik), cores = 8)
Q1_looic <- loo(Q1_log_lik, r_eff = Q1_r_eff)$estimates[3, 1:2]

Q1_looic
