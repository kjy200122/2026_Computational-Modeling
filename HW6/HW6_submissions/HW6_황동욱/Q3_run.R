
rm(list=ls())
setwd("~/Desktop/CCS Lab/Computational Modeling/숙제 2026/HW6/HW6_Dongwook")

library(rstan)

source("simulate_hw6_q3.R")

dat = read.table("simul_data_hw6_q3.txt", header=T, sep="\t")

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

output_q3 = stan("hw6_model1.stan", data = dataList,
                 pars = c("alpha", "beta"),
                 iter = 2000, warmup=1000, chains=2, cores=2)

print(output_q3)

# Extract posterior
parameters_q3 <- rstan::extract(output_q3)

alpha_mean_q3 = apply(parameters_q3$alpha, 2, mean)
alpha_sd_q3   = apply(parameters_q3$alpha, 2, sd)
beta_mean_q3  = apply(parameters_q3$beta, 2, mean)
beta_sd_q3    = apply(parameters_q3$beta, 2, sd)

# Scatter plots
par(mfrow=c(1,2))

plot(simul_pars$alpha, alpha_mean_q3,
     xlim=c(0, 0.5), ylim=c(0, 0.5),
     xlab="True alpha", ylab="Estimated alpha (posterior mean)",
     main="Q3: alpha recovery (200 subj, 100 trials)")
abline(0, 1, col="red")
arrows(x0=simul_pars$alpha, y0=alpha_mean_q3 - alpha_sd_q3,
       y1=alpha_mean_q3 + alpha_sd_q3,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean_q3,
     xlim=c(0, 4), ylim=c(0, 4),
     xlab="True beta", ylab="Estimated beta (posterior mean)",
     main="Q3: beta recovery (200 subj, 100 trials)")
abline(0, 1, col="red")
arrows(x0=simul_pars$beta, y0=beta_mean_q3 - beta_sd_q3,
       y1=beta_mean_q3 + beta_sd_q3,
       length=0.05, angle=90, code=3)

par(mfrow=c(1,1))
