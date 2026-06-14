rm(list=ls())
setwd("~/Desktop/CCS Lab/Computational Modeling/숙제 2026/HW6/HW6_Dongwook")

library(rstan)

source("simulate_hw6_model1.R")

dat = read.table("simul_data_hw6_model1.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)
N = length(allSubjs)
T = table(dat$subjID)[1]

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

# Run the NON-hierarchical Stan model
output_nh = stan("hw6_model1_nonhier.stan", data = dataList,
                 pars = c("alpha", "beta"),
                 iter = 2000, warmup=1000, chains=2, cores=2)

print(output_nh)

# Extract posterior
parameters_nh <- rstan::extract(output_nh)

alpha_mean_nh = apply(parameters_nh$alpha, 2, mean)
alpha_sd_nh   = apply(parameters_nh$alpha, 2, sd)
beta_mean_nh  = apply(parameters_nh$beta, 2, mean)
beta_sd_nh    = apply(parameters_nh$beta, 2, sd)

# Scatter plots
par(mfrow=c(1,2))

plot(simul_pars$alpha, alpha_mean_nh,
     xlim=c(0, 0.5), ylim=c(0, 0.5),
     xlab="True alpha", ylab="Estimated alpha (posterior mean)",
     main="Q2: alpha recovery (non-hierarchical)")
abline(0, 1, col="red")
arrows(x0=simul_pars$alpha, y0=alpha_mean_nh - alpha_sd_nh,
       y1=alpha_mean_nh + alpha_sd_nh,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean_nh,
     xlim=c(0, 4), ylim=c(0, 4),
     xlab="True beta", ylab="Estimated beta (posterior mean)",
     main="Q2: beta recovery (non-hierarchical)")
abline(0, 1, col="red")
arrows(x0=simul_pars$beta, y0=beta_mean_nh - beta_sd_nh,
       y1=beta_mean_nh + beta_sd_nh,
       length=0.05, angle=90, code=3)

par(mfrow=c(1,1))
